#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sun Feb  2 10:27:04 2020

@author: ligk2e
"""

import pandas as pd
import os
import task1_mod as mich  # current working directory should be concordant with task1_mod.py

def match_with_dbase(info_dict,df_ori,dbase):
    df_exam=[]
    df_back=[]
    for i in range(df_ori.shape[0]):
        exam,back=[],[]
        temp = mich.UID(df_ori,i)
        idn = list(temp.keys())[0]
        events_exam = list(temp.values())[0][0]
        events_back = list(temp.values())[0][1]
        for val in info_dict[idn][0][events_exam][1].values():
            for j in range(len(dbase)):
                if val in list(dbase.values())[j]:
                    exam.append(list(dbase.keys())[j])
                else:
                    pass
        for val in info_dict[idn][0][events_back][1].values():
            for j in range(len(dbase)):
                if val in list(dbase.values())[j]:
                    back.append(list(dbase.keys())[j])
                else:
                    pass
#        exam = query_longest(exam,dbase)
#        back = query_longest(back,dbase)        
        df_exam.append(exam)
        df_back.append(back)
    df_ori['exam_match']=df_exam
    df_ori['back_match']=df_back
    return df_ori
                

def query_longest(listq,dbase):
    if listq:
        max = 0
        transcriptID=''
        for item in listq:
            if dbase[item] > max:
                max = dbase[item]
                transcriptID = item
        listq=[dbase[item]]
    return listq  
    
        


if __name__ == '__main__':
#    df_whole = pd.read_csv('/Users/ligk2e/Desktop/project/PSI.AML__U2AF1-CV_vs_Healthy__U2AF1-CV.txt',
#                           sep='\t')
#    df_increase = df_whole[df_whole['dPSI']>0]
#    df_decrease = df_whole[df_whole['dPSI']<0]
#    df_increase.to_csv('/Users/ligk2e/Desktop/df_increase.txt',sep='\t',header=True,index=False)
#    df_decrease.to_csv('/Users/ligk2e/Desktop/df_decrease.txt',sep='\t',header=True,index=False)
    df_ori = pd.read_csv('df_decrease.txt',sep='\t')
    translation_de = mich.translate_junction('df_decrease.txt','\t')
    dbase = mich.make_dict('/project/SEQUENCE-protein-dbase_exoncomp.txt')
    result = match_with_dbase(translation_de,df_ori,dbase)
    result.to_csv('ban1.txt',sep='\t',header=True,index=False)
    
    