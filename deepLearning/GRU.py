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
    
    
# let's extract all possible contact site
seq1_contact = []
seq2_contact = []
for counter in dic1.values():
    for k,v in counter.items():
        if v > 10: seq1_contact.append(seq1[k])

for counter in dic2.values():
    for k,v in counter.items():
        if v > 10: seq2_contact.append(seq2[k])


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
import torch.nn as nn
import torch.nn.functional as F

from sklearn.metrics import roc_auc_score

#torch.manual_seed(0)

whole = []
for i in range(data.shape[0]):
    peptide = data.iloc[i]['peptide']
    cat = ''.join(seq1_contact) + peptide + ''.join(seq2_contact)
    whole.append(cat)
data_new = data.join(pd.Series(whole,name='whole'))
    

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
            onehot = torch.from_numpy(my_dataset.onehot_peptide(whole)) # length*20
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
            result[:,i] = my_dataset.onehot_aa(pep[i])
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





# let's build the model
class GRU_immuno(nn.Module):
    def __init__(self,seq_len,hidden_size,batch_size,input_size,num_layers,bi):
        super(GRU_immuno, self).__init__()
        self.seq_len = seq_len
        self.hidden_size = hidden_size
        self.batch_size = batch_size
        self.input_size = input_size
        self.num_layers = num_layers
        self.bi = bi
        self.gru3 = nn.GRU(input_size,hidden_size,num_layers=num_layers, bidirectional = bi,batch_first=True)
        self.attn = nn.Linear(2*hidden_size+1,1)   # bidirectional, so *2 , then add the next hidden size, in total, *3
        self.gru1 = nn.GRU(2*hidden_size,1,num_layers=1,bidirectional=False,batch_first=True)
        
    def forward(self,x):
        # print(x.size())  torch.Size([10, 191, 20])
        out,h_n = self.gru3(x.float(),self.init_hidden_state())
        next_init = torch.zeros(1*1,self.batch_size,1)
        weights= []
        for i in range(out.shape[1]):
            # print(torch.squeeze(out[:,i,:],dim=1).size()) torch.Size([10, 200])
            # print(torch.squeeze(next_init,dim=0).size()) torch.Size([10, 1])
            input_attn = torch.cat((torch.squeeze(out[:,i,:],dim=1),torch.squeeze(next_init,dim=0)),dim=1)
            # print(input_attn.size())  # torch.Size([10, 201])
            weights.append(self.attn(input_attn))
        normalized_weights = F.softmax(torch.cat(weights,dim=1),dim=1).unsqueeze(1)  # [batch_size,1,seq_len]
        context = torch.bmm(normalized_weights,out)     # [batch_size ,1 ,2*hidden_size]
        result,_ = self.gru1(context,next_init)   # [batch_size,1,1]
        # print(result.size())  torch.Size([10, 1, 1])
        # print(F.sigmoid(result.squeeze(1)).size())  torch.Size([10, 1])
        return torch.sigmoid(result.squeeze(1))
        
        
        
    def init_hidden_state(self):
        return torch.zeros((int(self.bi)+1)*self.num_layers,self.batch_size,self.hidden_size)
        
        
        
class Estimator(object):
    def __init__(self,model,optimizer,loss_f):
        self.model = model
        self.optimizer = optimizer
        self.loss_f = loss_f
        
    def training_one_epoch(self,dataLoader,batch_size):
        loss_list = []
        acc_list = []
        for idx,(X,y) in enumerate(dataLoader):
            X,y = X.float(),y.long()
            self.optimizer.zero_grad()
            y_pred = self.model(X)
            y_neg = torch.ones([batch_size,1]) - y_pred
            y_pred4crossEntropy = torch.cat([y_neg,y_pred],dim=1)  # [batch_size,2]
            # print(y_pred4crossEntropy.size(),y_pred.dtype)
            # print(y.size(),y.dtype)
            loss = self.loss_f(y_pred4crossEntropy,y)
            loss.backward()
            self.optimizer.step()
            loss_list.append(loss.item())
            #print(y_pred[:,0])
            #print((y_pred[:,0]>0.5).float())
            #print(y.float())

            accuracy = ((y_pred[:,0] > 0.5).float() == y.float()).sum().data.numpy()/batch_size
            #print(accuracy)
            acc_list.append(accuracy)
            
        return sum(loss_list)/len(loss_list),sum(acc_list)/len(acc_list)
    
    def fit(self,dataLoader,batch_size,epoch):
        self.model.train()
        for i in range(epoch):
            loss,acc = self.training_one_epoch(dataLoader,batch_size)
            print('Epoch {0}/{1} loss: {2:6.2f} - accuracy: {3:6.2f}'.format(i,epoch,loss,acc))
            
    def validation(self,dataset):  # for validation set
        val_loader = DataLoader(dataset,batch_size=len(dataset),shuffle=False)
        for idx,(X,y) in enumerate(val_loader):
            y_pred = self.model(X)
            y_pred_np = y_pred.data.numpy()
            auc = roc_auc_score(y_pred_np,y.data.numpy())
        return auc
    
    def prediction(self,dataset):  # for testing set
        test_loader = DataLoader(dataset,batch_size=10,shuffle=True,drop_last=True)
        for idx,(X,y) in enumerate(test_loader):
            y_pred = self.model(X)
            y_pred_np = (y_pred > 0.5).float().long().view(1,-1).squeeze(0).data.numpy()
            y_label = y.long().data.numpy()
            print(y_pred_np,y_label)
            auc = roc_auc_score(y_pred_np,y_label)
        return auc
                
            
#training_set, validation_set, testing_set = random_split(data_torch,(200,5,4))
training_set, testing_set = random_split(data_torch,(170,39))
        

# Starting to work
data_torch_training_loader = DataLoader(training_set,batch_size=10,shuffle=True,drop_last=True)
model = GRU_immuno(seq_len=191, hidden_size=100, batch_size=10, input_size=20, num_layers=3, bi=True)
clf = Estimator(model,optimizer=torch.optim.Adam(model.parameters(), lr=0.001, weight_decay = 0.001),loss_f=nn.CrossEntropyLoss())         
clf.fit(data_torch_training_loader,10,5)   
clf.prediction(testing_set)     
        
            
            
    




































