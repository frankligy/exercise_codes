#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Feb  3 10:15:26 2020

@author: ligk2e
"""
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

os.chdir('/Users/ligk2e/Desktop/IDA1/')

class decision_tree:
    def __init__(self,df,classinfo=None):
        self.df = df
        self.classinfo = classinfo
        
    def Split(self,test_size,random_state,y_colname):
        X = self.df.drop(y_colname,axis=1)
        Y = self.df[y_colname]
        X_train,X_test,Y_train,Y_test = train_test_split(X,y,
                            test_size=test_size,random_state=random_state)
        return X_train,X_test,Y_train,Y_test
    
    def Construction(X_train,Y_train,parameter=None):
        clf = tree.DecisionTreeClassifier(min_samples_leaf = parameter[])
        clf = clf.fit(X_train, Y_train)
    #    clf.predict([[2,2,2,2,2,2]]) # pay heed to the dimention, (1,6) ndarray
    #    clf.predict_proba([[2,2,2,2,2,2]])
        dot_data = StringIO(tree.export_graphviz(clf,out_file=None))
        graph = pydotplus.graph_from_dot_data(dot_data.getvalue())
        graph.write_png('tree{}.png'.format(parameter[]))
        return clf
    
    def Prediction(X_test,Y_test,clf):
        metrics_dict={}
        Y_pred = clf.predict(X_test)
        metrics_dict['accuracy'] = metrics.accuracy_score(Y_test,Y_pred))
        for group in classinfo.keys():
            metrics['precision_{}'.format(classinfo[group]) =metrics.precision_score(
                    Y_test,Y_pred,pos_label=classinfo[]))
            metrics['recall_{}'.format(classinfo[group]) =metrics.recall_score(
                    Y_test,Y_pred,pos_label=classinfo[]))
        return metrics_dict
   ''' 
    def draw_plot(metrics_dict):
        length = len(metrics_dict)
        style_set = {'-g','--r','-.k',':c','-m','--b','-y'}
        fig = plt.figure()
        ax = plt.axes()
        for key in metrics_dict.keys():
            ax.plot(parameter[],'-g',label='accuracy')
            ax.plot(leafnode_length,precision[class1],'--r',label='precision[class1]')
            ax.plot(leafnode_length,precision[class2],'-.k',label='precision[class2]')
            ax.plot(leafnode_length,recall[class1],':c',label='recall[class1]')
            ax.plot(leafnode_length,recall[class2],'-m',label='recall[class2]')
            ax.set(xlim=(0,55),ylim=(0,1),xlabel='min_records_leafnode',ylabel='metrics',title='metrics plot')
            ax.legend()
            fig.savefig(path)
            plt.close(fig)
      ''' 
                 
if __name__ == "__main__":
    df = pd.read_csv('Biomechanical_Data_2Classes.csv',sep=',')
    q1 = decision_tree(df,{'class1':'Abnormal','class2':'Normal'})
    X_train,X_test,Y_train,Y_test = q1.Split(80/230,1)
    accuracy,precision1,precision2,recall1,recall2 = [],[],[],[],[]
    for val in [5,15,25,40,50]: 
        clf = q1.Construction(X_train,Y_train,{'leafnode_record'= val})
        metrics_dict = q1.Prediction(X_test,Y_test,clf)
    
    
    
    
    
