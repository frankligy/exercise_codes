#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sat Feb  1 21:09:20 2020

@author: ligk2e
"""

# pip install scikit-learn

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


def split_train_test(train_num,seed):
    random.seed(seed)
    train = np.array(random.sample(range(310),train_num))
    test = np.delete(range(310),train)
    X_train = data2.values[train,:6]
    Y_train = data2.values[train,6]
    X_test = data2.values[test,:6]
    Y_test = data2.values[test,6]
    return X_train,X_test,Y_train,Y_test

# X2_train,X2_test,Y2_train,Y2_test = train_test_split(X,y,test_size=0.3,random_state=1)

def decision_tree2(class1,class2):
    accuracy=[]
    precision={class1:[],class2:[]}
    recall={class1:[],class2:[]}
    
    for i in range(5):
        clf = tree.DecisionTreeClassifier(min_samples_leaf = leafnode_length[i])
        clf = clf.fit(X2_train, Y2_train)
    #    clf.predict([[2,2,2,2,2,2]]) # pay heed to the dimention, (1,6) ndarray
    #    clf.predict_proba([[2,2,2,2,2,2]])
        dot_data = StringIO(tree.export_graphviz(clf,out_file=None))
        graph = pydotplus.graph_from_dot_data(dot_data.getvalue())
        graph.write_png('/Users/ligk2e/Desktop/IDA1/tree{}.png'.format(leafnode_length[i]))
        
        Y2_pred = clf.predict(X2_test)
        accuracy.append(metrics.accuracy_score(Y2_test,Y2_pred))
        precision[class1].append(metrics.precision_score(Y2_test,Y2_pred,
                 pos_label=class1))
        precision[class2].append(metrics.precision_score(Y2_test,Y2_pred,
                 pos_label=class2))
        recall[class1].append(metrics.recall_score(Y2_test,Y2_pred,
              pos_label=class1))
        recall[class2].append(metrics.recall_score(Y2_test,Y2_pred,
              pos_label=class2))
    return accuracy,precision,recall
    #    disp = metrics.plot_precision_recall_curve(clf,X2_test,Y2_test)
    #    disp.ax_.set_title('2-class Precisio-Recall curve with tree{}'.format(
    #            leafnode_length[i]))
        #try scikitplot.metrics.plot_precision_recall_curve could have more elegant PR curve
# draw the plot
def draw_dt2(class1,class2,accuracy,precision,recall,path):
    fig = plt.figure()
    ax = plt.axes()
    ax.plot(leafnode_length,accuracy,'-g',label='accuracy')
    ax.plot(leafnode_length,precision[class1],'--r',label='precision[class1]')
    ax.plot(leafnode_length,precision[class2],'-.k',label='precision[class2]')
    ax.plot(leafnode_length,recall[class1],':c',label='recall[class1]')
    ax.plot(leafnode_length,recall[class2],'-m',label='recall[class2]')
    ax.set(xlim=(0,55),ylim=(0,1),xlabel='min_records_leafnode',ylabel='metrics',title='metrics plot')
    ax.legend()
    fig.savefig(path)
    plt.close(fig)


if __name__ == '__main__':
    data2 = pd.read_csv('/Users/ligk2e/Desktop/Biomechanical_Data_2Classes.csv',sep=',')
    leafnode_length = [5,15,25,40,50]
    X2_train,X2_test,Y2_train,Y2_test = split_train_test(230,3)
    accuracy,precision,recall = decision_tree2('Abnormal','Normal')
    draw_dt2('Abnormal','Normal',accuracy,precision,recall,'/Users/ligk2e/Desktop/IDA1/metrics.svg')



