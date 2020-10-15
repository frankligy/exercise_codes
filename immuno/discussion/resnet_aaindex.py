'''
This script is going to use 34 HLA pseudo-sequence instead of paratope info,
still train the ResNet network, see if performance change

download all hla_prot.fasta from https://www.ebi.ac.uk/ipd/imgt/hla/download.html via FTP

1. deephlapan and netMHCpan pseudo-sequence they use, which is 34 of length
2. own improvement
'''

import pandas as pd
import numpy as np
import tensorflow as tf
import tensorflow.keras as keras
from tensorflow.keras import layers
import matplotlib.pyplot as plt
from sklearn.metrics import precision_recall_curve,roc_curve,auc,confusion_matrix


def chunk_clean_pseudo34():
    pseudo = pd.read_csv('immuno2/hla_deephla.txt', sep='\t', header=None)
    # delete all non-human alleles
    conds = []
    for i in range(pseudo.shape[0]):
        cond = True
        hla = pseudo[0].iloc[i]
        if not hla.startswith('HLA-'):
            cond = False
        conds.append(cond)
    pseudo_clean = pseudo.loc[conds]
    # change the format of hla
    clean_hla = [item[0:5] + '*' + item[5:7] + item[8:10] for item in pseudo_clean[0]]
    pseudo_clean[0] = clean_hla
    # test if length of pseudo sequence is 34 for all
    sum([1 if len(item) != 34 else 0 for item in pseudo_clean[1]])
    # write it to data folder for further use
    pseudo_clean.to_csv('immuno2/data/pseudo34_clean.txt', sep='\t', index=None)
    # I then change the column name to 'HLA' and 'pseudo'

def chuck_valid_training():
    ori = pd.read_csv('immuno2/data/ineo_testing_filter910_new.txt', sep='\t')
    hla = pd.read_csv('immuno2/data/pseudo34_clean.txt',sep='\t')
    # make sure all training dataset has valid hla allele
    pseudo = set(hla['HLA'].tolist())
    conds = []
    for i in range(ori.shape[0]):
        cond = True
        col1 = ori['HLA'].iloc[i]
        if not col1 in pseudo:
            cond = False
        conds.append(cond)
    ori = ori.loc[conds]
    ori.to_csv('immuno2/data/ineo_testing_filter910_valid.txt',sep='\t',index=None)

def aaindex(peptide,after_pca):

    amino = 'ARNDCQEGHILKMFPSTWYV-'
    matrix = np.transpose(after_pca)   # [12,21]
    encoded = np.empty([len(peptide), 12])  # (seq_len,12)
    for i in range(len(peptide)):
        query = peptide[i]
        if query == 'X': query = '-'
        query = query.upper()
        encoded[i, :] = matrix[:, amino.index(query)]

    return encoded

def peptide_data_aaindex(peptide,after_pca):   # return numpy array [10,12,1]
    length = len(peptide)
    if length == 10:
        encode = aaindex(peptide,after_pca)
    elif length == 9:
        peptide = peptide[:5] + '-' + peptide[5:]
        encode = aaindex(peptide,after_pca)
    encode = encode.reshape(encode.shape[0], encode.shape[1], -1)
    return encode


def dict_inventory(inventory):
    dicA, dicB, dicC = {}, {}, {}
    dic = {'A': dicA, 'B': dicB, 'C': dicC}

    for hla in inventory:
        type_ = hla[4]  # A,B,C
        first2 = hla[6:8]  # 01
        last2 = hla[8:]  # 01
        try:
            dic[type_][first2].append(last2)
        except KeyError:
            dic[type_][first2] = []
            dic[type_][first2].append(last2)

    return dic


def rescue_unknown_hla(hla, dic_inventory):
    type_ = hla[4]
    first2 = hla[6:8]
    last2 = hla[8:]
    big_category = dic_inventory[type_]
    #print(hla)
    if not big_category.get(first2) == None:
        small_category = big_category.get(first2)
        distance = [abs(int(last2) - int(i)) for i in small_category]
        optimal = min(zip(small_category, distance), key=lambda x: x[1])[0]
        return 'HLA-' + str(type_) + '*' + str(first2) + str(optimal)
    else:
        small_category = list(big_category.keys())
        distance = [abs(int(first2) - int(i)) for i in small_category]
        optimal = min(zip(small_category, distance), key=lambda x: x[1])[0]
        return 'HLA-' + str(type_) + '*' + str(optimal) + str(big_category[optimal][0])






def hla_data_aaindex(hla_dic,hla_type,after_pca):    # return numpy array [34,12,1]
    try:
        seq = hla_dic[hla_type]
    except KeyError:
        hla_type = rescue_unknown_hla(hla_type,dic_inventory)
        seq = hla_dic[hla_type]
    encode = aaindex(seq,after_pca)
    encode = encode.reshape(encode.shape[0], encode.shape[1], -1)
    return encode

def construct_aaindex(ori,hla_dic,after_pca):
    series = []
    for i in range(ori.shape[0]):
        peptide = ori['peptide'].iloc[i]
        hla_type = ori['HLA'].iloc[i]
        immuno = np.array(ori['immunogenecity'].iloc[i]).reshape(1,-1)   # [1,1]

        encode_pep = peptide_data_aaindex(peptide,after_pca)    # [10,12]

        encode_hla = hla_data_aaindex(hla_dic,hla_type,after_pca)   # [46,12]
        series.append((encode_pep, encode_hla, immuno))
    return series

def hla_df_to_dic(hla):
    dic = {}
    for i in range(hla.shape[0]):
        col1 = hla['HLA'].iloc[i]  # HLA allele
        col2 = hla['pseudo'].iloc[i]  # pseudo sequence
        dic[col1] = col2
    return dic

def pull_peptide_aaindex(dataset):
    result = np.empty([len(dataset),10,12,1])
    for i in range(len(dataset)):
        result[i,:,:,:] = dataset[i][0]
    return result


def pull_hla_aaindex(dataset):
    result = np.empty([len(dataset),46,12,1])
    for i in range(len(dataset)):
        result[i,:,:,:] = dataset[i][1]
    return result


def pull_label_aaindex(dataset):
    result = np.empty([len(dataset),1])
    for i in range(len(dataset)):
        result[i,:] = dataset[i][2]
    return result

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


class CNN_peptide_aaindex(layers.Layer):
    def __init__(self):
        super(CNN_peptide_aaindex,self).__init__()
        self.conv = layers.Conv2D(filters=16,kernel_size=(3,12),strides=(1,1))
        self.block1 = ResBlock(16,(2,1))
        self.block2 = ResBlock(32,(2,1))
        self.block3 = ResBlock(64,(2,1))

    def call(self,x):    # (10,21,1)
        out = self.conv(x)   # (8,1,16)
        out = self.block1(out)   # (4,1,32)
        out = self.block2(out)   # (2,1,64)
        out = self.block3(out)   # (1,1,128)
        return out


class CNN_MHC_aaindex(layers.Layer):
    def __init__(self):
        super(CNN_MHC_aaindex,self).__init__()
        self.conv = layers.Conv2D(filters=16,kernel_size=(15,12),strides=(1,1)) # (32,1,16)
        self.block1 = ResBlock(16, (2, 1))    # (16,1,32)
        self.block2 = ResBlock(32, (2, 1))    # (8,1,64)
        self.block3 = ResBlock(64, (2, 1))    # (4,1,128)
        self.conv_add = layers.Conv2D(filters=128,kernel_size=(4,1),strides=(1,1))
        self.bn = layers.BatchNormalization()


    def call(self, x):
        out = self.conv(x)
        #print(out.shape)
        out = self.block1(out)
        #print(out.shape)
        out = self.block2(out)
        out = self.block3(out)
        out = keras.activations.relu(self.bn(self.conv_add(out)))   # (1,1,128)
        return out


class model_aaindex(keras.Model):
    def __init__(self):
        super(model_aaindex,self).__init__()
        self.br_pep = CNN_peptide_aaindex()
        self.br_mhc = CNN_MHC_aaindex()
        self.flatten = layers.Flatten()
        self.fc1 = layers.Dense(128,activation='relu')
        self.fc2 = layers.Dense(1,activation='sigmoid')

    def call(self,input):
        x1,x2 = input[0],input[1]  # x1: (10,12,1)    x2: (46,12,1)
        out1 = self.flatten(self.br_pep(x1))
        out2 = self.flatten(self.br_mhc(x2))
        out = layers.concatenate([out1,out2])
        out = self.fc1(out)
        out = self.fc2(out)
        return out

    def model(self):
        x1 = keras.Input(shape=(10,12,1))
        x2 = keras.Input(shape=(46,12,1))
        return keras.Model(inputs=[x1,x2],outputs=self.call([x1,x2]))


def draw_ROC(y_true,y_pred):

    fpr,tpr,_ = roc_curve(y_true,y_pred,pos_label=1)
    area_mine = auc(fpr,tpr)
    fig = plt.figure()
    lw = 2
    plt.plot(fpr, tpr, color='darkorange',
            lw=lw, label='ROC curve (area = %0.2f)' % area_mine)
    plt.plot([0, 1], [0, 1], color='navy', lw=lw, linestyle='--')
    plt.xlim([0.0, 1.0])
    plt.ylim([0.0, 1.05])
    plt.xlabel('False Positive Rate')
    plt.ylabel('True Positive Rate')
    plt.title('Receiver operating characteristic example')
    plt.legend(loc="lower right")
    plt.show()

def draw_PR(y_true,y_pred):
    precision,recall,_ = precision_recall_curve(y_true,y_pred,pos_label=1)
    area_PR = auc(recall,precision)
    baseline = np.sum(np.array(y_true) == 1) / len(y_true)

    plt.figure()
    lw = 2
    plt.plot(recall,precision, color='darkorange',
            lw=lw, label='PR curve (area = %0.2f)' % area_PR)
    plt.plot([0, 1], [baseline, baseline], color='navy', lw=lw, linestyle='--')
    plt.xlim([0.0, 1.0])
    plt.ylim([0.0, 1.05])
    plt.xlabel('Recall')
    plt.ylabel('Precision')
    plt.title('PR curve example')
    plt.legend(loc="lower right")
    plt.show()

def draw_history(history):
    plt.subplot(211)
    plt.title('Loss')
    plt.plot(history.history['loss'], label='train')
    plt.plot(history.history['val_loss'], label='validation')
    plt.legend()
    # plot accuracy during training
    plt.subplot(212)
    plt.title('Accuracy')
    plt.plot(history.history['accuracy'], label='train')
    plt.plot(history.history['val_accuracy'], label='validation')
    plt.legend()
    plt.show()

if __name__ == '__main__':
    after_pca = np.loadtxt('immuno2/data/after_pca.txt')
    ori = pd.read_csv('immuno2/data/shuffle_all_filter910_valid.txt', sep='\t')
    hla = pd.read_csv('immuno2/data/hla2paratopeTable_aligned.txt',sep='\t')
    hla_dic = hla_df_to_dic(hla)
    inventory = list(hla_dic.keys())
    dic_inventory = dict_inventory(inventory)

    dataset = construct_aaindex(ori,hla_dic,after_pca)
    input1 = pull_peptide_aaindex(dataset)
    input2 = pull_hla_aaindex(dataset)
    label = pull_label_aaindex(dataset)

    resnet_aaindex = model_aaindex()
    resnet_aaindex.compile(
        loss='binary_crossentropy',
        optimizer=keras.optimizers.Adam(lr=0.0001),
        metrics=['accuracy'])

    # let's do a simple train, validation split
    array = np.arange(len(dataset))
    train_index = np.random.choice(array,int(len(dataset)*0.9),replace=False)
    valid_index = [item for item in array if item not in train_index]

    input1_train = input1[train_index]
    input1_valid = input1[valid_index]
    input2_train = input2[train_index]
    input2_valid = input2[valid_index]
    label_train = label[train_index]
    label_valid = label[valid_index]

    callback_val = keras.callbacks.EarlyStopping(monitor='val_loss', patience=20,restore_best_weights=False)
    callback_train = keras.callbacks.EarlyStopping(monitor='loss',patience=2,restore_best_weights=False)
    history = resnet_aaindex.fit(
        x=[input1_train,input2_train],   # feed a list into
        y=label_train,
        validation_data = ([input1_valid,input2_valid],label_valid),
        batch_size=128,
        epochs=200,
        class_weight = {0:0.2,1:0.8},   # I have 20% positive and 80% negative in my training data
        callbacks = [callback_val,callback_train])

    # testing-virus dataset
    ori_test = pd.read_csv('immuno2/data/ineo_testing_filter910_valid.txt',sep='\t')
    dataset_test = construct_aaindex(ori_test,hla_dic,after_pca)
    input1_test = pull_peptide_aaindex(dataset_test)
    input2_test = pull_hla_aaindex(dataset_test)
    label_test = pull_label_aaindex(dataset_test)

    pred = resnet_aaindex.predict(x=[input1_test,input2_test])
    draw_ROC(label_test,pred)
    draw_PR(label_test,pred)

    # neoantigen_mannual dataset
    ori_test = pd.read_csv('immuno2/data/mannual_cancer_testing_fiter910.txt',sep='\t')
    dataset_test = construct_aaindex(ori_test,hla_dic,after_pca)
    input1_test = pull_peptide_aaindex(dataset_test)
    input2_test = pull_hla_aaindex(dataset_test)
    label_test = pull_label_aaindex(dataset_test)

    pred = resnet_aaindex.predict(x=[input1_test,input2_test])
    draw_ROC(label_test,pred)
    draw_PR(label_test,pred)


    ## Machine learning algorithm






























