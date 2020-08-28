import matplotlib.pyplot as plt
import tensorflow as tf
from tensorflow import keras
from tensorflow.keras import layers
from sklearn.metrics import roc_curve, auc, precision_recall_curve, confusion_matrix, f1_score,accuracy_score
import pandas as pd
import numpy as np
from seperateCNN import *
from utils import *
import shelve

def draw_combined_ROC(arrayP,arrayT):
    fig = plt.figure()
    colormap = ['red','green','blue','orange']
    legend = ['iedb','logistic regression','deephlapan','deepimmuno']
    for i in range(len(arrayP)):
        fpr,tpr,_ = roc_curve(arrayT[i],arrayP[i],pos_label=1)
        area = auc(fpr, tpr)
        lw = 2
        plt.plot(fpr, tpr, color=colormap[i],
                 lw=lw, label='{0} (area = {1:0.2f})'.format(legend[i],area))
    plt.plot([0, 1], [0, 1], color='navy', lw=lw, linestyle='--')
    plt.xlim([0.0, 1.0])
    plt.ylim([0.0, 1.05])
    plt.xlabel('False Positive Rate')
    plt.ylabel('True Positive Rate')
    plt.title('Receiver operating characteristic example')
    plt.legend(loc="lower right")
    plt.show()

def draw_combined_PR(arrayP,arrayT):
    fig = plt.figure()
    colormap = ['red','green','blue','orange']
    legend = ['iedb','logistic regression','deephlapan','deepimmuno']
    for i in range(len(arrayP)):
        precision,recall,_ = precision_recall_curve(arrayT[i],arrayP[i],pos_label=1)
        area = auc(recall, precision)
        lw = 2
        plt.plot(recall, precision, color=colormap[i],
                 lw=lw, label='{0} (area = {1:0.2f})'.format(legend[i],area))
    plt.plot([0, 1], [0.12, 0.12], color='navy', lw=lw, linestyle='--')
    plt.xlim([0.0, 1.0])
    plt.ylim([0.0, 1.05])
    plt.xlabel('Recall')
    plt.ylabel('Precision')
    plt.title('PR curve example')
    plt.legend(loc="upper right")
    plt.show()

if __name__ == '__main__':
    # load testing dataset
    ori_test = pd.read_csv('data/ineo_testing_filter910_new.txt',sep='\t')  # shuffle_validation_filter910.txt # ineo_testing_filter910_new.txt
    hla = pd.read_csv('hla2paratopeTable_aligned.txt',sep='\t',header=None,names=['hla','paratope'])
    inventory = hla['hla']
    dic_inventory = dict_inventory(inventory)
    dataset_test = construct(ori_test, hla, dic_inventory)
    input1_test = pull_peptide(dataset_test)
    input2_test = pull_hla(dataset_test)
    label_test = pull_label(dataset_test)




    # result from deepimmuno
    seperateCNNmodel = seperateCNN()
    seperateCNNmodel.load_weights('secondFilter32_epoch150/')
    result = seperateCNNmodel.predict(x=[input1_test,input2_test])

    # result from baseline logistic regression
    s = shelve.open('logistic')
    y_pred = s['y_pred']
    Y_test = s['Y_test']
    s.close()
    y_pred = y_pred[:,1]
    Y_test = Y_test.astype(np.int32)

    # result from deephlapan
    df1 = pd.read_csv('deephlapan/ineo_testing_new.txt', sep='\t')
    df2 = pd.read_csv('deephlapan/ineo_testing_new_final_predicted_result.csv')
    y_deephlapan = df1['immunogenecity'].values
    y_pred_deephlapan = df2['immunogenic score'].values

    # result from IEDB
    iedb = pd.read_csv('IEDB.csv')
    y_iedb = iedb['label'].tolist()
    y_pred_iedb = iedb['score'].tolist()

    # let's draw the figure
    arrayP = [y_pred_iedb,y_pred,y_pred_deephlapan,result[:,1]]
    arrayT = [y_iedb,Y_test,y_deephlapan,label_test]
    draw_combined_ROC(arrayP,arrayT)
    draw_combined_PR(arrayP,arrayT)



