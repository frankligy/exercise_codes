#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Jul 15 18:07:09 2020

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
    data_clean2 = data_clean.loc[data_clean[cond]]
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
    data_clean2.update({'HLA':[new]})
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
        
    
class my_dataset(Dataset):
    def __init__(self,original):
        self.original = original
        self.new_data = self.padded_to_max_len()   # it is a torch
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
            onehot = torch.from_numpy(my_dataset.encode_peptide(whole)) # length*20
            #print(onehot,onehot.size(),onehot[0,:])
            tmp.append(onehot)
        #print(tmp[0][0,:])
        stacked_equal_length = pad_sequence(tmp,batch_first=True)  #[number of item,length,20]
        
        return stacked_equal_length
    
    def get_Y(self):
        return torch.from_numpy(self.original['immunogenecity'].values).view(1,-1).float()
            
        
    
    @staticmethod
    def encode_peptide(pep):
        result = np.zeros([20,len(pep)])
        for i in range(result.shape[1]):
            
            result[:,i] = my_dataset.onehot_aa(pep[i])

        return np.transpose(result)
    
    @staticmethod
    def onehot_aa(a):
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
        mat = torch.eye(len(transform))
        result = mat[transform[a],:]
        
        return result
        
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
        for i in matrix.shape[0]:
            for j in matrix.shape[1]:
                matrix[i,j] = dic[(transform_rev[i],transform_rev[j])] 
                
                
        return matrix[:,transform[a]]
        
        
        
        
        

        



class Encoder(nn.Module):
    def __init__(self,input_size,hidden_size,num_layers):
        super(Encoder,self).__init__()
        self.input_size = input_size
        self.hidden_size = hidden_size
        self.num_layers = num_layers
        
        self.gru = nn.GRU(self.input_size,self.hidden_size,num_layers=self.num_layers,
                          bidirectional=True,batch_first=True)
        
    def forward(self,x): # shape of x: [batch,seq_len,input_size]

        out,hn = self.gru(x)
        return out,hn
    
    

class Decoder(nn.Module):
    def __init__(self,hidden_size,num_layers=1):
        super(Decoder,self).__init__()
        
        # self.encoder_out = encoder_out
        # self.encoder_hn = encoder_hn
        
        self.output_size = hidden_size   # encoder_hidden_size
        self.input_size = hidden_size*2 +  self.output_size  # 2 * encoder_hidden_size
        self.hidden_size = self.output_size
        
        self.gru = nn.GRU(self.input_size,self.hidden_size,num_layers=1,bidirectional=False,batch_first=True)
        
        self.energy = nn.Linear(self.input_size,1)
        self.Softmax = nn.Softmax(dim=2)
        self.relu = nn.ReLU()
        
    def forward(self,x,encoder_out,encoder_hn): # shape of x : [batch_size,1,self.output_size]
    
        weights = []   # store softmax(Eij) value
        #print('encoder_out:',encoder_out.size())
        for i in range(encoder_out.shape[1]):
            tmp = x.transpose(0,1)   # change the shape of x to [1,batch_size,output_size] for torch.cat
            # print('tmp:',tmp.size())
            # print('encoder_out[:,i,:]:',encoder_out[:,i,:].size())
            
            # then remember the shape of encoder_out[i] is: [1,batch_size,2*encoder_hidden_size]
            E = self.relu(self.energy(torch.cat((tmp,encoder_out[:,i,:].unsqueeze(0)),dim=2)))  # [1,batch_size,1]
            weights.append(E)
            
        stacked_weights = torch.cat(weights,dim=2)  # [1,batch_size,seq_len]
        normalized_weights = self.Softmax(stacked_weights) # didn't change dimensions
        context_vector = torch.bmm(normalized_weights.transpose(0,1),encoder_out)  
        # [batch,1,seq_Len] * [batch,seq_len,2*encoder_hidden_size] = [batch,1,2*encoder_hidden_size]
        
        out,hn = self.gru(torch.cat((x,context_vector),dim=2))
        
        return out,hn
    
class GRU_immuno(nn.Module):
    def __init__(self,encoder,decoder,target_len,hidden_size,linear_hidden_size,p=0.25):
        super(GRU_immuno,self).__init__()
        self.encoder = encoder
        self.decoder = decoder
        self.target_len = target_len
        self.hidden_size = hidden_size
        self.linear_hidden_size = linear_hidden_size
        
        self.fc1 = nn.Linear(self.target_len*self.hidden_size,linear_hidden_size)
        self.dropout = nn.Dropout(p)
        self.relu = nn.ReLU()
        self.fc2 = nn.Linear(self.linear_hidden_size,2)
        self.Sigmoid = nn.Sigmoid()
        
    def forward(self,input_data,target0):  # [batch,seq_len,input_size]  # target0 : [batch,1,outsize/hiddensize]
         
        encoder_out,encoder_hn = self.encoder(input_data)
        #print('encoder_out:',encoder_out.size())
        
        running_hn = encoder_hn
        x = target0
        
        batch_size = encoder_out.shape[0]
        hidden_size = encoder_hn.shape[2]
        
        result = torch.zeros([self.target_len,batch_size,hidden_size])
        for i in range(self.target_len):
            decoder_out,decoder_hn = self.decoder(x,encoder_out,running_hn)
            # dimention of decoder_out: [batch,1,hidden]
            # print('decoder',decoder_out.size())
            # print('result',result[i].size())
            result[i] = decoder_out.transpose(0,1).squeeze(0)
            x = decoder_out
            running_hn = decoder_hn
        
        # shape of result : [target_len,batch_size,hidden_size]
            
        result_go = result.reshape(result.shape[1],-1)
        out = self.relu(self.dropout(self.fc1(result_go)))
        out = self.fc2(out)
        out = self.Sigmoid(out)   # out: [batch,2]
            
        return out
            
            
        
        
        
            

                
            
def count_parameters(model):
    table = PrettyTable(['Modules','Parameters'])
    total_params = 0
    for name,parameter in model.named_parameters():
        if not parameter.requires_grad: continue
        param = parameter.numel()
        table.add_row([name,param])
        total_params += param
    print(table)
    print('Total Trainable Params:{0}'.format(total_params))

    

def shuffleTensor(t):
    idx = torch.randperm(t.nelement())
    t = t.view(-1)[idx].view(t.size())
    return t

def permuteAlongAxis(t,axis):
    idx = torch.randperm(t.size()[axis])
    t = t[idx,:,:]
    return t
    
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
        
    
        

    
    loader = zip(chunks_X_list,chunks_y_list)
    return loader
        
        
        
    
    
    

    
    
    
    
    
    
    


if __name__ == '__main__':
    
    start_time = time.time()
    
    data = pd.read_excel('/Users/ligk2e/Desktop/NeoAntigenWorkflow/immunogenecity/data/data.xlsx')
    data_clean = clean_data_frame(data)
    
    hla_seq = read_hla('/Users/ligk2e/Desktop/NeoAntigenWorkflow/immunogenecity/hla_seq/hla_seq.txt')
    inventory = pd.read_csv('/Users/ligk2e/Desktop/NeoAntigenWorkflow/immunogenecity/imgt_contact/inventory.txt',sep='\n',header=None)[0].tolist() 
    dic_inventory = dict_inventory(inventory)   
    data_whole = get_whole_seq('/Users/ligk2e/Desktop/NeoAntigenWorkflow/immunogenecity/imgt_contact',data_clean)
    
    
    data_torch = my_dataset(data_whole)
    #training_set, validation_set, testing_set = random_split(data_torch,(200,5,4))
    training_set, testing_set = random_split(data_torch,(10240,1057))
            
    device = torch.device('cuda' if torch.cuda.is_available() else 'cpu')
    
    # Starting to work
    # data_torch_training_loader = DataLoader(training_set,batch_size=512,shuffle=True,drop_last=True)
    # data_torch_testing_loader = DataLoader(testing_set,batch_size=512,shuffle=True,drop_last=True)
    
    data_torch_training_loader = balancedBinaryLoader(training_set,batch_size=64)
    
    
    encoder = Encoder(input_size=20,hidden_size=50,num_layers=3).to(device)
    decoder = Decoder(hidden_size=50,num_layers=1).to(device)
    
    model = GRU_immuno(encoder, decoder, target_len=30, hidden_size=50, linear_hidden_size=50).to(device)
    
    optimizer = torch.optim.Adam(model.parameters(), lr=0.0001)
    scheduler = scheduler = torch.optim.lr_scheduler.ReduceLROnPlateau(
        optimizer, factor=0.1, patience=2, verbose=True)
    # if it observe a non-decreasing loss, give you another (patience-1) more chances, if still not decrease, will reduce 
    # learning rate to factor*learning rate
    loss_f=nn.CrossEntropyLoss()
    target0 = torch.randn([64,1,50]).to(device)
    
    num_epochs = 5
    for epoch in range(num_epochs):

        for i in data_torch_training_loader:
            loss_list = []
            acc_list = []
            X,y = i[0].float().to(device),i[1].long().to(device)
            #print(X,y)
            optimizer.zero_grad()
            
            y_pred = model(X,target0)
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
        


    
    end_time = time.time()
    print('Consumed {0} seconds'.format(end_time-start_time))
    
    

    
    
    

    
    
    
    
    
    
    