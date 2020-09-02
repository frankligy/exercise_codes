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
        self.maxpool = layers.MaxPool2D(pool_size=pool_size,strides=(1,1))



    def call(self,x):
        out = keras.activations.relu(self.bn1(self.conv1(x)))
        out = keras.activations.relu(self.bn2(self.conv2(out)))
        out = keras.activations.relu(self.bn3(self.conv3(out)))
        identity_map = self.identity(x)
        out = out + identity_map
        out = self.maxpool(out)

        return out


class CNN_peptide(layers.Layer):
    def __init__(self):
        super(CNN_peptide,self).__init__()
        self.conv = layers.Conv2D(filters=16,kernel_size=(3,21),strides=(1,1))
        self.block1 = ResBlock(16,(2,1))
        self.block2 = ResBlock(32,(2,1))
        self.block3 = ResBlock(64,(2,1))

    def call(self,x):
        out = self.conv(x)
        out = self.block1(out)
        out = self.block2(out)
        out = self.block3(out)
        return out


class CNN_MHC(layers.Layer):
    def __init__(self):
        super(CNN_MHC,self).__init__()
        self.conv = layers.Conv2D(filters=16,kernel_size=(15,21),strides=(1,1))
        self.block1 = ResBlock(16, (9, 1))
        self.block2 = ResBlock(32, (9, 1))
        self.block3 = ResBlock(64, (9, 1))

    def call(self, x):
        out = self.conv(x)
        out = self.block1(out)
        out = self.block2(out)
        out = self.block3(out)
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
        x1,x2 = input[0],input[1]
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

    ResLikeCNN = model()


    ResLikeCNN.compile(
        loss=keras.losses.SparseCategoricalCrossentropy(from_logits=True),
        optimizer=keras.optimizers.Adam(lr=0.001),
        metrics=['accuracy']
    )
    ResLikeCNN.fit(
        x=[input1,input2],
        y=label,
        batch_size=512,
        epochs=50,
        #class_weight = {0:0.2,1:0.8}
    )







