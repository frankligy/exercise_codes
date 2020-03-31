#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sun Mar 29 16:48:52 2020

@author: ligk2e
"""

import os
os.chdir("/Users/ligk2e/Desktop/IDA2/q1/")
import numpy as np
import pandas as pd
from sklearn.model_selection import train_test_split
from sklearn import svm
from sklearn import metrics
import matplotlib.pyplot as plt
#import matpotlib.colors as colors
#import matplotlib.cm as cmx
import xlrd


class Data():
    
    def __init__(self,data,label):
        self.data = data
        self.label = label

    def getData(self):
        return self.data
    
    def getLabel(self):
        return self.label
    
    def setData(self,newData):
        self.data = newData
        
    def setLabel(self,newLabel):
        self.label = newLabel
    
    def svm(self,gamma=0.1,C=1.0):
        X_train,X_test,Y_train,Y_test = train_test_split(self.data,self.label,test_size = 0.25,random_state = 109)
        clf = svm.SVC(kernel='rbf',C=C,gamma=gamma)
        clf.fit(X_train,Y_train)
        Y_pred = clf.predict(X_test)
        confusionMatrix = metrics.confusion_matrix(Y_test,Y_pred)
        precision = metrics.precision_score(Y_test,Y_pred)
        recall = metrics.recall_score(Y_test,Y_pred)
        accuracy = metrics.accuracy_score(Y_test,Y_pred)
        print('Radial-basis kernel: gamma = {0}, C={1}:\n'.format(gamma,C))
        print('Confusion matrix is:\n{0}'.format(confusionMatrix))
        print('Precision is:{}\n'.format(precision))
        print('Recall is:{}\n'.format(recall))
        print('Accuracy is:{}\n'.format(accuracy))
        return clf,precision, recall, accuracy
        
        
def plotMetrics(paramArray,precisionArray,recallArray,accuracyArray,path):
    fig = plt.figure()
    ax = plt.axes()
    ax.plot(paramArray,precisionArray,'*g',label='precision')
    ax.plot(paramArray,recallArray,'1k',label='recall')
    ax.plot(paramArray,accuracyArray,'sm',label='accuray')
    ax.axvline(0.4,label='highest accuray: C=0.4')
    ax.set(xlim=(0,2),ylim=(0,1),xlabel='value of regularization parameter C', ylabel = 'Value of metrics',
           title='metrics plot')
    ax.legend()
    plt.show()
    fig.savefig(path)
    plt.close(fig)     

def predOnGrid(clf,path):
    gridData = np.zeros((10000,2),dtype=np.int16)
    for i in range(100):  #row
        for j in range(100):  #column
            gridData[100*i+j] = [i+1,j+1]
    print(gridData,gridData.shape)
    pred = clf.predict(gridData)
    # scatter plot, color coding by label method 1, using pandas groupby method
    df=pd.DataFrame({'x':gridData[:,0],'y':gridData[:,1],'label':pred})   
    groups = df.groupby('label')
    '''
    label x y
    1    4  5
    1    3  7
    0    4  0
    ...
    0    5  6
    '''
    fig = plt.figure()
    ax = plt.axes()
    for name,group in groups:
        ax.plot(group.x,group.y,marker='o',linestyle='',label=name)
    ax.legend()
    ax.set_title('Classification on 100*100 grid')
    plt.show()
    fig.savefig(path)
    plt.close(fig)
    return df
    
    
    
            


if __name__ == "__main__":
    data = pd.read_excel('HW2Data.xlsx',sheet_name='Sheet1',names=['x','y','blank','label'])
    data = data.drop(['blank'],axis=1)
    hwData = Data(data[['x','y']].values,data['label'].values)
    precisionArray,recallArray,accuracyArray,paramArray = [],[],[],[]
    for C in [0.1,0.2,0.3,0.4,0.5,0.6,0.7,1,1.5]:
        clf, precision, recall, accuracy = hwData.svm(0.1,C)
        precisionArray.append(precision)
        recallArray.append(recall)
        accuracyArray.append(accuracy)
        paramArray.append(C)
    plotMetrics(paramArray,precisionArray,recallArray,accuracyArray,'svmMetrics.pdf')   
    
    
    clf,precision,recall,accuray = hwData.svm(0.1,0.4)
    
    
    df = predOnGrid(clf,'gridPlot.pdf')
    

# Python painting: https://www.learnpyqt.com/courses/custom-widgets/bitmap-graphics/
# import pyqt5    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
        