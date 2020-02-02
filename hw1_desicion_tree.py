#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sat Feb  1 21:09:20 2020

@author: ligk2e
"""

# pip install scikit-learn

from sklearn import tree
from sklearn.model_selection import train_test_split
from sklearn import metrics
from sklearn.externals.six import StringIO
import pydotplus
import os
import numpy as np
import pandas as pd
import random
import matplotlib.pyplot as plt
import graphviz

data2 = pd.read_csv('/Users/ligk2e/Desktop/Biomechanical_Data_2Classes.csv',sep=',')
data3 = pd.read_csv('/Users/ligk2e/Desktop/Biomechanical_Data_3Classes.csv',sep=',')
random.seed(3)
train = np.array(random.sample(range(310),230))
test = np.delete(range(310),train)
X2_train = data2.values[train,:6]
Y2_train = data2.values[train,6]
X2_test = data2.values[test,:6]
Y2_test = data2.values[test,6]

# X2_train,X2_test,Y2_train,Y2_test = train_test_split(X,y,test_size=0.3,random_state=1)


accuracy=[]
precision={'Abnormal':[],'Normal':[]}
recall={'Abnormal':[],'Normal':[]}
leafnode_length = [5,15,25,40,50]
for i in range(5):
    clf = tree.DecisionTreeClassifier(min_samples_leaf = leafnode_length[i])
    clf = clf.fit(X2_train, Y2_train)
#    clf.predict([[2,2,2,2,2,2]]) # pay heed to the dimention, (1,6) ndarray
#    clf.predict_proba([[2,2,2,2,2,2]])
    
#    a = tree.plot_tree(clf.fit(X2_train,Y2_train))
#    
#    dot_data = tree.export_graphviz(clf,out_file=None)
#    graph_list.append(graphviz.Source(dot_data))
    
    #dot_data = StringIO()
    dot_data = StringIO(tree.export_graphviz(clf,out_file=None))
    graph = pydotplus.graph_from_dot_data(dot_data.getvalue())
    graph.write_png('/Users/ligk2e/Desktop/IDA1/tree{}.png'.format(leafnode_length[i]))
# all the graph has been saved to /Users/ligk2e/Desktop/IDA
    
    Y2_pred = clf.predict(X2_test)
    accuracy.append(metrics.accuracy_score(Y2_test,Y2_pred))
    precision['Abnormal'].append(metrics.precision_score(Y2_test,Y2_pred,
             pos_label="Abnormal"))
    precision['Normal'].append(metrics.precision_score(Y2_test,Y2_pred,
             pos_label="Normal"))
    recall['Abnormal'].append(metrics.recall_score(Y2_test,Y2_pred,
          pos_label="Abnormal"))
    recall['Normal'].append(metrics.recall_score(Y2_test,Y2_pred,
          pos_label="Normal"))
    disp = metrics.plot_precision_recall_curve(clf,X2_test,Y2_test)
    disp.ax_.set_title('2-class Precisio-Recall curve with tree{}'.format(
            leafnode_length[i]))
    #try scikitplot.metrics.plot_precision_recall_curve could have more elegant PR curve








