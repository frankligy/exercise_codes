import os
import tensorflow as tf
from tensorflow import keras
from tensorflow.keras import layers,regularizers
import pandas as pd
from utils import *
from sklearn.metrics import roc_curve, auc, precision_recall_curve, confusion_matrix, f1_score,accuracy_score
import matplotlib.pyplot as plt
import numpy as np

from seperateCNN import *


class ResBlock(layers.Layer):
    def __init__(self,in_channel,pool_size):
        super(ResBlock,self).__init__()
        intermediate_channel = in_channel
        out_channel = in_channel * 2
        self.conv1 = layers.Conv2D(filters=intermediate_channel,kernel_size=(1,1),strides=(1,1),padding='same')
        self.bn1 = layers.BatchNormalization()
        self.conv2 = layers.Conv2D(filters=intermediate_channel,kernel_size=(3,1),strides=(1,1),padding='same')
        self.bn2 = layers.BatchNormalization()
        self.conv3 = layers.Conv2D(filters=out_channel,kernel_size=(1,1),strides=(1,1),padding='same')
        self.bn3 = layers.BatchNormalization()
        self.identity = layers.Conv2D(filters=out_channel,kernel_size=(1,1),strides=(1,1))
        self.maxpool = layers.MaxPool2D(pool_size=pool_size,strides=pool_size)



    def call(self,x):
        out = keras.activations.relu(self.bn1(self.conv1(x)))   # (8,1,16)
        out = keras.activations.relu(self.bn2(self.conv2(out)))  # (8,1,16)
        out = keras.activations.relu(self.bn3(self.conv3(out)))   # (8,1,32)
        identity_map = self.identity(x)   # change original input (8,1,16)  --> (8,1,32)
        out = out + identity_map    # (8,1,32)
        out = self.maxpool(out)    # (4,1,32)

        return out


class CNN_peptide(layers.Layer):
    def __init__(self):
        super(CNN_peptide,self).__init__()
        self.conv = layers.Conv2D(filters=16,kernel_size=(3,21),strides=(1,1))
        self.block1 = ResBlock(16,(2,1))
        self.block2 = ResBlock(32,(2,1))
        self.block3 = ResBlock(64,(2,1))

    def call(self,x):    # (10,21,1)
        out = self.conv(x)   # (8,1,16)
        out = self.block1(out)   # (4,1,32)
        out = self.block2(out)   # (2,1,64)
        out = self.block3(out)   # (1,1,128)
        return out


class CNN_MHC(layers.Layer):
    def __init__(self):
        super(CNN_MHC,self).__init__()
        self.conv = layers.Conv2D(filters=16,kernel_size=(15,21),strides=(1,1)) # (32,1,16)
        self.block1 = ResBlock(16, (2, 1))    # (16,1,32)
        self.block2 = ResBlock(32, (2, 1))    # (8,1,64)
        self.block3 = ResBlock(64, (2, 1))    # (4,1,128)
        self.conv_add = layers.Conv2D(filters=128,kernel_size=(4,1),strides=(1,1))
        self.bn = layers.BatchNormalization()


    def call(self, x):
        out = self.conv(x)
        out = self.block1(out)
        out = self.block2(out)
        out = self.block3(out)
        out = keras.activations.relu(self.bn(self.conv_add(out)))   # (1,1,128)
        return out


class model(keras.Model):
    def __init__(self):
        super(model,self).__init__()
        self.br_pep = CNN_peptide()
        self.br_mhc = CNN_MHC()
        self.flatten = layers.Flatten()
        self.fc1 = layers.Dense(128,activation='relu')
        self.fc2 = layers.Dense(2,activation='softmax')

    def call(self,input):
        x1,x2 = input[0],input[1]  # x1: (10,21,1)    x2: (46,21,1)
        out1 = self.flatten(self.br_pep(x1))
        out2 = self.flatten(self.br_mhc(x2))
        out = layers.concatenate([out1,out2])
        out = self.fc1(out)
        out = self.fc2(out)
        return out




if __name__ == '__main__':
    ori = pd.read_csv('data/shuffle_training_test.txt',sep='\t')  # shuffle_validation_filter910.txt # ineo_testing_filter910.txt
    hla = pd.read_csv('hla2paratopeTable_aligned.txt',sep='\t',header=None,names=['hla','paratope'])
    inventory = hla['hla']
    dic_inventory = dict_inventory(inventory)

    dataset = construct(ori, hla, dic_inventory)   # [ (10,21,1),(46,21,1),(1,1)   ]
    input1 = pull_peptide(dataset)
    input2 = pull_hla(dataset)
    label = pull_label(dataset)

    '''
    input1: ndarray (28581,10,21,1)
    input2: ndarray (28581,46,21,1)
    label : ndarray (28581,1)
    '''

    ResLikeCNN = model()


    ResLikeCNN.compile(
        loss=keras.losses.SparseCategoricalCrossentropy(from_logits=True),
        optimizer=keras.optimizers.Adam(lr=0.001),
        metrics=['accuracy']
    )
    ResLikeCNN.fit(
        x=[input1,input2],   # feed a list into
        y=label,
        batch_size=512,
        epochs=20,
        class_weight = {0:0.4,1:0.6}   # I have 20% positive and 80% negative in my training data
    )


# now let's test in external dataset
    ori_test = pd.read_csv('data/ineo_testing_filter910_new.txt',sep='\t')  # shuffle_validation_filter910.txt # ineo_testing_filter910_new.txt
    dataset_test = construct(ori_test, hla, dic_inventory)
    input1_test = pull_peptide(dataset_test)
    input2_test = pull_hla(dataset_test)
    label_test = pull_label(dataset_test)
    #seperateCNNmodel.evaluate(x=[input1_test,input2_test],y=label_test,batch_size=512)
    result = ResLikeCNN.predict(x=[input1_test,input2_test])
    hard = [1 if i > 0.5 else 0 for i in result[:,1]]
    confusion_matrix(label_test,hard)
    f1_score(label_test,hard)
    accuracy_score(label_test,hard)
    draw_ROC(label_test,result[:,1])
    draw_PR(label_test,result[:,1])




