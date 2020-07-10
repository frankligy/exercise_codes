#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sun Mar 29 16:48:52 2020

@author: ligk2e
"""

import os
os.chdir("/Users/ligk2e/Desktop/IDA2/q2/")
import numpy as np
import pandas as pd
from sklearn.model_selection import train_test_split
from sklearn import svm
from sklearn import metrics
import matplotlib.pyplot as plt
from statistics import mean,stdev
from math import sqrt
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
    
    def kMeans(self,k=3,n=1,mode='q1'):
        from sklearn.cluster import KMeans
        if mode == 'q1':
            km = KMeans(n_clusters = k, init = 'random', n_init = n, max_iter = 300, tol = 1e-04) # km will belong to KMeans class
        elif mode == 'q2':
            km = KMeans(n_clusters = k, init = 'random', n_init = n, max_iter = 300, tol = 1e-04,random_state=0)
        km.fit_predict(self.data)   # label of each point, array(n,1), store in km.labels_
        return km.cluster_centers_, km.inertia_,km.labels_   # trailing single underscore is to avoid confliting with built-in variables
    
    def kmeansPlot(self, labels, centroids,path):
        palette = ['green','pink','orange','blue','black'] # green, pink, red, blue, black # hexadecimal color
        labelsUnique = list(set(labels))  #{0,1,2,3,4} set object is a little bit tricky, it is unsubscriptable
        clusterSSE = []
        for i in range(len(labelsUnique)):            
            sum_ = clusterSpecificSSE(self.data[labels == labelsUnique[i],:],centroids[i,:])
            clusterSSE.append(sum_)
            plt.scatter(self.data[labels == labelsUnique[i],0],self.data[labels==labelsUnique[i],1],s = 50,
                        c = palette[i],marker = 'o', edgecolor = 'black', label = 'cluster{0}, SSE{0} = {1}'.format(labelsUnique[i],clusterSSE[i])) # all the points
        plt.scatter(centroids[:,0],centroids[:,1],s =250, marker = '*', c ='red', edgecolor = 'black', label = 'centroids')
        plt.legend(bbox_to_anchor=(1.04,1))  # nice! Remember that!
        plt.grid()
        plt.savefig(path,bbox_inches='tight') # prior to plt.show() since plt.show() close the image
        plt.show()
     
    def intuitiveClustering(self):
        self.intuitiveLabel = []
        for i in range(self.data.shape[0]):
            x,y = self.data[i,0],self.data[i,1]
            if x >= 0 and x < 15 and y >=0 and y<60: self.intuitiveLabel.append(3)
            elif x >= 15 and x < 55 and y >=0 and y < 60: self.intuitiveLabel.append(0)
            elif x >= 55 and x < 100 and y >=0 and y < 60: self.intuitiveLabel.append(2)
            elif x >=0 and x <= 100 and y >= 60 and y < 80: self.intuitiveLabel.append(1)
            elif x >= 0 and x <= 100 and y >= 80 and y <=100: self.intuitiveLabel.append(4)
            
def clusterSpecificSSE(array2d,list_):  #(n,2)
    sum_ = 0
    for i in range(array2d.shape[0]):
        L2Norm = sqrt((array2d[i,0] - list_[0])**2 + (array2d[i,1] - list_[1])**2)  # L2 norm is Euclidean distance and L1 norm is Manhattan distance
        sum_ += L2Norm
    return Round2Precision(sum_,1)
        
    
             
def Round2Precision(value,precision:int=0,mode:str=''): # default argument with specified type
    assert precision >= 0 # if true, continue, otherwise raise assertError, using for self-check
    value *= 10 ** precision # if you wanna round by precision, have to do that
    method = round   # round will base on >.5 or <.5
    if mode.lower() == 'up': 
        method = math.ceil     # always round up
    elif mode.lower() == 'down':
        method = math.floor   # always round down
    answer = '{0:.{1}f}'.format(method(value)/10**precision,precision)   
    return float(answer)         

def meanPlusPrecision(list):
    return Round2Precision(mean(list),1) 

def sdPlusPrecision(list):
    return Round2Precision(stdev(list),1)      
        
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
    data = pd.read_excel('HW2Data(1).xlsx',sheet_name='Sheet1',names=['x','y','blank','label'])
    data = data.drop(['blank'],axis=1)
    hwData = Data(data[['x','y']].values,data['label'].values)
###### svm portion
#    precisionArray,recallArray,accuracyArray,paramArray = [],[],[],[]
#    for C in [0.1,0.2,0.3,0.4,0.5,0.6,0.7,1,1.5]:
#        clf, precision, recall, accuracy = hwData.svm(0.1,C)
#        precisionArray.append(precision)
#        recallArray.append(recall)
#        accuracyArray.append(accuracy)
#        paramArray.append(C)
#    plotMetrics(paramArray,precisionArray,recallArray,accuracyArray,'svmMetrics.pdf')   
#    
#    
#    clf,precision,recall,accuray = hwData.svm(0.1,0.4)
#    
#    
#    df = predOnGrid(clf,'gridPlot.pdf')
    

# Python painting: https://www.learnpyqt.com/courses/custom-widgets/bitmap-graphics/
# import pyqt5    
    
    
    
##### Kmeans portion
    sseArrayList = []
    for k in [3,5,7,9,11]:
        sseArray = []
        for j in range(6):
            _,sse,_ = hwData.kMeans(k)
            sseArray.append(Round2Precision(sse,1))
        sseArrayList.append(sseArray)

    averageSSE = list(map(meanPlusPrecision,sseArrayList))
    sdSSE = list(map(sdPlusPrecision,sseArrayList))
    maxSSE = list(map(max,sseArrayList))
    minSSE = list(map(min,sseArrayList))
    print('Total SSE for each k value in each run are shown below: \n')
    print(sseArrayList,'\n')
    print('---------------------------------------------')
    print('Average SSE for each k value are shown below:\n')
    print(averageSSE,'\n')
    print('---------------------------------------------')
    print('Standard Deviation for each k value are shown below:\n')
    print(sdSSE,'\n')
    print('---------------------------------------------')
    print('Max SSE for each k value are shown below:\n')
    print(maxSSE,'\n')
    print('---------------------------------------------')
    print('Min SSE for each k value are shown below:\n')
    print(minSSE,'\n')
    
    
    centroids, sse5, labels = hwData.kMeans(5,6,'q2')

    hwData.kmeansPlot(labels,centroids,'kmeansPlot.pdf')
    
    from sklearn.metrics import adjusted_rand_score
    hwData.intuitiveClustering()
    RI = adjusted_rand_score(hwData.intuitiveLabel,labels)
    print('rand index is {0}\n'.format(RI))
    
    
    
    
    
    
    
    
    
    
    
        