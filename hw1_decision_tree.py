from sklearn import tree
import os
import numpy as np
import pandas as pd
import random

data2 = pd.read_csv('C:\\Users\\ligk2E\\Desktop\\coding\\DT\\Biomechanical_Data_2Classes.csv',sep=',')
data3 = pd.read_csv('C:\\Users\\ligk2E\\Desktop\\coding\\DT\\Biomechanical_Data_3Classes.csv',sep=',')
train = np.array(random.sample(range(310),230))
test = np.delete(range(310),train)
X2_train = data2.values[train-1,:6]
Y2_train = data2.values[train-1,6]
X2_test = data2.values[test-1,:6]
Y2_test = data2.values[test-1,6]

clf = tree.DecisionTreeClassifier(min_samples_leaf = 5)
clf = clf.fit(X2_train, Y2_train)