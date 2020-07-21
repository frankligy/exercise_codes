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

from torch.utils.data import Dataset, DataLoader, random_split
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
            onehot = torch.from_numpy(my_dataset.onehot_peptide(whole)) # length*20
            #print(onehot,onehot.size(),onehot[0,:])
            tmp.append(onehot)
        #print(tmp[0][0,:])
        stacked_equal_length = pad_sequence(tmp,batch_first=True)  #[number of item,length,20]
        
        return stacked_equal_length
    
    def get_Y(self):
        return torch.from_numpy(self.original['immunogenecity'].values).view(1,-1).float()
            
        
    
    @staticmethod
    def onehot_peptide(pep):
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
            
            y_pred = self.model(X,torch.randn([batch_size,1,hidden_size]).to(device))
            loss = self.loss_f(y_pred,y)
            loss.backward()
            self.optimizer.step()
            loss_list.append(loss.item())
            
            num_correct = 0
            num_samples = 0
            _,predictions = y_pred.max(1)

            num_correct += (predictions == y).sum()  # will generate a 0-D tensor, tensor(49), float() to convert it

            num_samples  += predictions.size(0)

            acc_list.append(float(num_correct)/float(num_samples)*100)
            
        return sum(loss_list)/len(loss_list),sum(acc_list)/len(acc_list)
    
    def fit(self,dataLoader,batch_size,epoch,hidden_size):
        self.model.train()
        for i in range(epoch):
            loss,acc = self.training_one_epoch(dataLoader,batch_size,hidden_size)
            self.scheduler.step(loss)
            print('Epoch {0}/{1} loss: {2:6.2f} - accuracy{3:6.2f}%'.format(i+1,epoch,loss,acc))
            
            
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
    data_torch_training_loader = DataLoader(training_set,batch_size=64,shuffle=True,drop_last=True)
    data_torch_testing_loader = DataLoader(testing_set,batch_size=64,shuffle=True,drop_last=True)
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
    
    num_epochs = 10
    for epoch in range(num_epochs):

        for idx,(X,y) in enumerate(data_torch_training_loader):
            loss_list = []
            acc_list = []
            X,y = X.float().to(device),y.long().to(device)
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
    
    

    
    
    

    
    
    
    
    
    
    