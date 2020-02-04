#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Feb  4 10:41:54 2020

@author: ligk2e
"""

from sklearn import tree
#from sklearn.model_selection import train_test_split
from sklearn import metrics
from sklearn.externals.six import StringIO
import pydotplus
import os
import numpy as np
import pandas as pd
import random
import matplotlib.pyplot as plt
import graphviz

def split_train_test(train_num,seed,df):
    random.seed(seed)
    train = np.array(random.sample(range(310),train_num))
    test = np.delete(range(310),train)
    X_train = df.values[train,:6]
    Y_train = df.values[train,6]
    X_test = df.values[test,:6]
    Y_test = df.values[test,6]
    return X_train,X_test,Y_train,Y_test

def decision_tree3(class1,class2,class3):
    accuracy=[]
    precision=[]
    recall=[]
    
    for i in range(5):
        clf = tree.DecisionTreeClassifier(min_samples_leaf = leafnode_length[i])
        clf = clf.fit(X2_train, Y2_train)
    #    clf.predict([[2,2,2,2,2,2]]) # pay heed to the dimention, (1,6) ndarray
    #    clf.predict_proba([[2,2,2,2,2,2]])
        dot_data = StringIO(tree.export_graphviz(clf,out_file=None))
        graph = pydotplus.graph_from_dot_data(dot_data.getvalue())
        graph.write_png('/Users/ligk2e/Desktop/IDA1/q2/tree{}.png'.format(leafnode_length[i]))
        
        Y2_pred = clf.predict(X2_test)
        accuracy.append(metrics.accuracy_score(Y2_test,Y2_pred))
        precision.append(metrics.precision_score(Y2_test,Y2_pred,average=None))
        recall.append(metrics.recall_score(Y2_test,Y2_pred,average=None))
    precision=extract(precision,'Hernia','Spondylolisthesis','Normal')
    recall=extract(recall,'Hernia','Spondylolisthesis','Normal')
    return accuracy,precision,recall

def extract(listq,class1,class2,class3):
    class1_a,class2_a,class3_a = [],[],[]
    for item in listq:
        class1_a.append(item[0])
        class2_a.append(item[1])
        class3_a.append(item[2])
    dict = {}
    dict[class1] = class1_a
    dict[class2] = class2_a
    dict[class3] = class3_a
    return dict

def draw_dt3(class1,class2,class3,accuracy,precision,recall,path):
    fig = plt.figure()
    ax = plt.axes()
    ax.plot(leafnode_length,accuracy,'-g',label='accuracy')
    ax.plot(leafnode_length,precision[class1],'--r',label='precision[{}]'.format(class1))
    ax.plot(leafnode_length,precision[class2],'-.k',label='precision[{}]'.format(class2))
    ax.plot(leafnode_length,precision[class3],'-y',label='precision[{}]'.format(class3))
    ax.plot(leafnode_length,recall[class1],':c',label='recall[{}]'.format(class1))
    ax.plot(leafnode_length,recall[class2],'-m',label='recall[{}]'.format(class2))
    ax.plot(leafnode_length,recall[class3],'-b',label='recall[{}]'.format(class3))
    ax.set(xlim=(0,55),ylim=(0,1),xlabel='min_records_leafnode',ylabel='metrics',title='metrics plot')
    ax.legend(fontsize=7)
    fig.savefig(path)
    plt.close(fig)
    
if __name__ == '__main__':
    data3 = pd.read_csv('/Users/ligk2e/Desktop/IDA1/Biomechanical_Data_3Classes.csv',sep=',')
    leafnode_length = [5,15,25,40,50]
    X2_train,X2_test,Y2_train,Y2_test = split_train_test(230,3,data3)
    accuracy,precision,recall = decision_tree3('Hernia','Spondylolisthesis','Normal')
    draw_dt3('Hernia','Spondylolisthesis','Normal',accuracy,precision,recall,'/Users/ligk2e/Desktop/IDA1/q2/metrics.svg')
    
    
    
    
    