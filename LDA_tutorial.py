#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Jun 29 22:39:44 2020

@author: ligk2e
"""


'''
python syntax:
1. np.set_printoptions(precision=4)
2. sns.set()  
3. wine = load_wine() 
4. df = X.join(pd.Series(y, name='class'))  # add a column
5. for c, rows in df.groupby('class'): 
        for index, row in rows.iterrows():  
6.  class_feature_means[c] = rows.mean()   # df.mean() return a Series, mean for each column
7. n = len(df.loc[df['class'] == c].index)   # filter by a column
8. pairs = [(np.abs(eigen_values[i]), eigen_vectors[:,i]) for i in range(len(eigen_values))] # tie them together
    pairs = sorted(pairs, key=lambda x: x[0], reverse=True)     
9. le = LabelEncoder()   # good to know
    y = le.fit_transform(df['class'])

.........
https://towardsdatascience.com/linear-discriminant-analysis-in-python-76b8b17817c2

premises:
1. normal distribution    
2. equal variance for each class

'''

from sklearn.datasets import load_wine
import pandas as pd
import numpy as np
np.set_printoptions(precision=4)
from matplotlib import pyplot as plt
import seaborn as sns
sns.set()       # very useful
from sklearn.preprocessing import LabelEncoder
from sklearn.tree import DecisionTreeClassifier
from sklearn.model_selection import train_test_split
from sklearn.metrics import confusion_matrix

wine = load_wine()   # wine is an object, haven't been a dataframe

X = pd.DataFrame(wine.data, columns=wine.feature_names)   
y = pd.Categorical.from_codes(wine.target, wine.target_names)
df = X.join(pd.Series(y, name='class'))  # add a column


class_feature_means = pd.DataFrame(columns=wine.target_names)
for c, rows in df.groupby('class'):   # c: 'class0', str   # rows: a sub-df, DataFrame                                     
    class_feature_means[c] = rows.mean()   # df.mean() return a Series, mean for each column
    
within_class_scatter_matrix = np.zeros((13,13))    
for c, rows in df.groupby('class'):
    rows = rows.drop(['class'], axis=1)
    s = np.zeros((13,13))
    for index, row in rows.iterrows():   # index is index, row is a Series object for each row
        x, mc = row.values.reshape(13,1), class_feature_means[c].values.reshape(13,1)  # to column vector
        s += (x - mc).dot((x - mc).T)  # numpy.outer(x,y)
        within_class_scatter_matrix += s
        


feature_means = df.mean()
between_class_scatter_matrix = np.zeros((13,13))
for c in class_feature_means: # we can iterate a df's column like this
    n = len(df.loc[df['class'] == c].index)   # filter by a column
    mc, m = class_feature_means[c].values.reshape(13,1), feature_means.values.reshape(13,1)
    between_class_scatter_matrix += n * (mc - m).dot((mc - m).T)
    
eigen_values, eigen_vectors = np.linalg.eig(np.linalg.inv(within_class_scatter_matrix).dot(between_class_scatter_matrix))
    
pairs = [(np.abs(eigen_values[i]), eigen_vectors[:,i]) for i in range(len(eigen_values))] # tie them together
pairs = sorted(pairs, key=lambda x: x[0], reverse=True)
for pair in pairs:
    print(pair[0])
    
eigen_value_sums = sum(eigen_values)
print('Explained Variance')
for i, pair in enumerate(pairs):
    print('Eigenvector {}: {}'.format(i, (pair[0]/eigen_value_sums).real))
    
w_matrix = np.hstack((pairs[0][1].reshape(13,1), pairs[1][1].reshape(13,1))).real  # combine (a,b) horizontally, if a and b are column vector, it looks like cbind                       
X_lda = np.array(X.dot(w_matrix))    # matrix dot

le = LabelEncoder()   # good to know
y = le.fit_transform(df['class'])

plt.xlabel('LD1')
plt.ylabel('LD2')
plt.scatter(
    X_lda[:,0],
    X_lda[:,1],
    c=y,
    cmap='rainbow',
    alpha=0.7,
    edgecolors='b'
)


from sklearn.discriminant_analysis import LinearDiscriminantAnalysis
lda = LinearDiscriminantAnalysis()
X_lda = lda.fit_transform(X, y)

lda.explained_variance_ratio_   # get leading eigenvalues

lda.predict_proba([[-3.6880e-01,  1.1966e+00,  2.1059e+00,  4.3845e-01, -1.0087e-02,
         2.6221e+00, -7.9606e+00, -9.0429e+00,  9.5294e-02,  1.9352e+00,
        -5.9296e+00, -4.9254e+00, -7.1364e-03]])
    
    
    
    
    
    
    
    
    
    
    
    
    
