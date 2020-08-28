import os
import tensorflow as tf
from tensorflow import keras
from tensorflow.keras import layers,regularizers
import pandas as pd
from utils import *
from sklearn.metrics import roc_curve, auc, precision_recall_curve, confusion_matrix, f1_score,accuracy_score
import matplotlib.pyplot as plt
import numpy as np


# construct the model use functional API
def seperateCNN():
    input1 = keras.Input(shape=(10, 21, 1))
    input2 = keras.Input(shape=(46, 21, 1))

    x = layers.Conv2D(filters=16, kernel_size=(2, 21))(input1)
    x = layers.BatchNormalization()(x)
    x = keras.activations.relu(x)
    x = layers.MaxPool2D(pool_size=(2, 1), strides=(1, 1))(x)
    x = layers.Conv2D(filters=32, kernel_size=(2, 1))(x)
    x = layers.BatchNormalization()(x)
    x = keras.activations.relu(x)
    x = layers.MaxPool2D(pool_size=(2, 1), strides=(1, 1))(x)
    x = layers.Flatten()(x)
    x = keras.Model(inputs=input1, outputs=x)

    y = layers.Conv2D(filters=16, kernel_size=(15, 21))(input2)
    y = layers.BatchNormalization()(y)
    y = keras.activations.relu(y)
    y = layers.MaxPool2D(pool_size=(5, 1), strides=(1, 1))(y)
    y = layers.Conv2D(filters=32,kernel_size=(15,1))(y)
    y = layers.BatchNormalization()(y)
    y = keras.activations.relu(y)
    y = layers.MaxPool2D(pool_size=(5, 1),strides=(1,1))(y)
    y = layers.Flatten()(y)
    y = keras.Model(inputs=input2,outputs=y)

    combined = layers.concatenate([x.output,y.output])
    # output_bias0 = keras.initializers.Constant(np.log([(28581-5496)/5496]))
    # output_bias1 = keras.initializers.Constant(np.log([5496/(28581-5496)]))
    z = layers.Dense(128,activation='relu')(combined)
    #z = layers.Dropout(0.3)(z)
    z = layers.Dense(2,activation='softmax')(z)

    model = keras.Model(inputs=[input1,input2],outputs=z)
    return model



def pull_peptide(dataset):
    result = np.empty([len(dataset),10,21,1])
    for i in range(len(dataset)):
        result[i,:,:,:] = dataset[i][0]
    return result


def pull_hla(dataset):
    result = np.empty([len(dataset),46,21,1])
    for i in range(len(dataset)):
        result[i,:,:,:] = dataset[i][1]
    return result


def pull_label(dataset):
    result = np.empty([len(dataset),1])
    for i in range(len(dataset)):
        result[i,:] = dataset[i][2]
    return result









if __name__ == '__main__':
    ori = pd.read_csv('data/shuffle_training_test.txt',sep='\t')  # shuffle_validation_filter910.txt # ineo_testing_filter910.txt
    hla = pd.read_csv('hla2paratopeTable_aligned.txt',sep='\t',header=None,names=['hla','paratope'])
    inventory = hla['hla']
    dic_inventory = dict_inventory(inventory)

    dataset = construct(ori, hla, dic_inventory)   # [ (10,21,1),(46,21,1),(1,1)   ]
    input1 = pull_peptide(dataset)
    input2 = pull_hla(dataset)
    label = pull_label(dataset)
    seperateCNNmodel = seperateCNN()
    seperateCNNmodel.load_weights('secondFilter32_epoch150/')
    seperateCNNmodel.compile(
        loss=keras.losses.SparseCategoricalCrossentropy(),
        optimizer=keras.optimizers.Adam(lr=0.0001),
        metrics=['accuracy']
    )
    seperateCNNmodel.fit(
        x=[input1,input2],
        y=label,
        batch_size=512,
        epochs=200,
        class_weight = {0:0.2,1:0.8}
    )

    # now let's test in external dataset
    ori_test = pd.read_csv('data/ineo_testing_filter910_new.txt',sep='\t')  # shuffle_validation_filter910.txt # ineo_testing_filter910_new.txt
    dataset_test = construct(ori_test, hla, dic_inventory)
    input1_test = pull_peptide(dataset_test)
    input2_test = pull_hla(dataset_test)
    label_test = pull_label(dataset_test)
    seperateCNNmodel.evaluate(x=[input1_test,input2_test],y=label_test,batch_size=512)
    result = seperateCNNmodel.predict(x=[input1_test,input2_test])
    hard = [1 if i > 0.5 else 0 for i in result[:,1]]
    confusion_matrix(label_test,hard)
    f1_score(label_test,hard)
    accuracy_score(label_test,hard)
    draw_ROC(label_test,result[:,1])
    draw_PR(label_test,result[:,1])

    # let's save the current model
    seperateCNNmodel.save_weights('secondFilter32_epoch150/')






