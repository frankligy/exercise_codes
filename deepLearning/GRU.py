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
        
    def forward(self,input_data,target0):  # [batch,seq_len,input_size]
         
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
            
            
        
        
        
            
            
        



        
        
        
class Estimator(object):
    def __init__(self,model,optimizer,scheduler,loss_f,device):
        self.model = model
        self.optimizer = optimizer
        self.scheduler = scheduler
        self.loss_f = loss_f
        self.device = device
        
    def training_one_epoch(self,dataLoader,batch_size,hidden_size):
        loss_list = []
        acc_list = []
        for idx,(X,y) in enumerate(dataLoader):
            X,y = X.float().to(device),y.long().to(device)
            self.optimizer.zero_grad()
            y_pred = self.model(X,torch.randn([batch_size,1,hidden_size]))
            loss = self.loss_f(y_pred,y)
            loss.backward()
            self.optimizer.step()
            loss_list.append(loss.item())
            
            num_correct = 0
            num_samples = 0
            _,predictions = y_pred.max(1)
            num_correct += (predictions == y).sum()
            num_samples  += predictions.size(0)
            acc_list.append(num_correct/num_samples)
            
        return sum(loss_list)/len(loss_list),sum(acc_list)/len(acc_list)
    
    def fit(self,dataLoader,batch_size,epoch,hidden_size):
        self.model.train()
        for i in range(epoch):
            loss,acc = self.training_one_epoch(dataLoader,batch_size,hidden_size)
            self.scheduler.step(loss)
            print('Epoch {0}/{1} loss: {2:6.2f} - accuracy{3:6.2f}'.format(i+1,epoch,loss,acc))
            
            
    def accuracy(self,loader):
        num_correct = 0
        num_samples = 0
        self.model.eval()
        
        with torch.no_grad():
            for x, y in loader:
                x = x.to(device=device)
                y = y.to(device=device)
    
                scores = self.model(x)
                _, predictions = scores.max(1)  # first is value, second is index, we need index
                num_correct += (predictions == y).sum()
                num_samples += predictions.size(0)
            print('Got {0}/{1} with accuracy {0}/{1}*100'.format(num_correct,num_samples))
        self.model.train()
            
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
        
device = torch.device('cuda' if torch.cuda.is_available() else 'cpu')

# Starting to work
data_torch_training_loader = DataLoader(training_set,batch_size=64,shuffle=True,drop_last=True)
data_torch_testing_loader = DataLoader(testing_set,batch_size=64,shuffle=True,drop_last=True)
encoder = Encoder(input_size=20,hidden_size=200,num_layers=3).to(device)
decoder = Decoder(hidden_size=200,num_layers=1).to(device)

model = GRU_immuno(encoder, decoder, target_len=50, hidden_size=200, linear_hidden_size=200).to(device)

optimizer = torch.optim.Adam(model.parameters(), lr=0.001, weight_decay = 0.001)
scheduler = scheduler = torch.optim.lr_scheduler.ReduceLROnPlateau(
    optimizer, factor=0.1, patience=2, verbose=True)
# if it observe a non-decreasing loss, give you another (patience-1) more chances, if still not decrease, will reduce 
# learning rate to factor*learning rate

clf = Estimator(model,optimizer=optimizer,
                scheduler = scheduler,
                loss_f=nn.CrossEntropyLoss(),device=device) 
        
clf.fit(data_torch_training_loader,batch_size=64,epoch=5,hidden_size=200)   
   
        
            
            
    




































