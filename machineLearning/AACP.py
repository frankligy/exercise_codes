#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sun Jun 28 17:05:49 2020

@author: ligk2e
"""
import pandas as pd
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

seq = seq1 + seq2 + seq3

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

# machine learning

# split training and testing set
X_train,X_test,Y_train,Y_test = train_test_split(X,Y,test_size =0.2,random_state=255)

# Random Forest
clf = RandomForestClassifier(min_samples_split=5,random_state =255)
clf.fit(X_train,Y_train)
clf.predict(X_test)
clf.predict_proba(X_test)
clf.predict_log_proba(X_test)
clf.score(X_test,Y_test)   # accurary: 0.57
clf.get_params()

scores = cross_val_score(clf,X,Y,cv=5)   # cross-validation, it's decent
'''
array([0.66666667, 0.61904762, 0.64285714, 0.76190476, 0.63414634])
'''


# ROC
fpr,tpr,thresholds = metrics.roc_curve(Y_test,clf.predict(X_test),pos_label=1)
metrics.auc(fpr,tpr)
metrics.roc_auc_score(Y_test,clf.predict(X_test))

# PR curve
precision,recall,thresholds = metrics.precision_recall_curve(Y_test,clf.predict(X_test))
metrics.plot_precision_recall_curve(clf,X_test,Y_test)

# Matthews correlation
corr = metrics.matthews_corrcoef(Y_test,clf.predict(X_test))

# SVM
clf = SVC(probability=True,random_state=1)
clf.fit(X_train,Y_train)
clf.predict(X_test)

scores = cross_val_score(clf,X,Y,cv=5)

# Neural Network
with open('/Users/ligk2e/Desktop/immunogenecity/X.p','wb') as f:
    pickle.dump(X,f)
    
with open('/Users/ligk2e/Desktop/immunogenecity/Y.p','wb') as f:
    pickle.dump(Y,f)
    
        
        
        
        
    
    
    









































        
    


    