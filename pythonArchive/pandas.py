#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Apr 28 10:31:00 2020

@author: ligk2e
"""

import pandas as pd
import numpy as np
import matplotlib.pyplot as plt

S1 = pd.Series(['a','b','c'],index=[1,3,4])  # Series object
S1.index  # return index, it is a Index object
S1.values   # return value, it is a narray object

Df = pd.DataFrame([['a','A'],['b','B'],['c','C']],columns=['lower','upper'],index=[0,1,2])
Df = pd.DataFrame({'lower':['a','b','c'],'upper':['A','B','C'],'wow':[1,2,3]},index=[0,1,2])

Df.index   # return index, it is a Index object
Df.columns   # return index, it is a Index object
Df.head(2)   
Df.shape
Df.info()    # could give information column-wise, useful for inspect NULL value
Df.describe()  # obtain some summary column-wise
Df.isnull()
Df.dropna(axis=1,how='all')  # eelete column if this column are all NULL, please look at the documentation
Df.fillna(0)
Df.drop_duplicates()

Df['lower'].dtype
Df['lower'].astype('float')

Df.index = [1,2,3]  # rename the index
Df.rename(index={1:'I',2:'II'})   # specify how to rename the index
Df.rename(columns={'lower':'lower1'})  # specify how to rename the columns

Df = Df[['lower','upper']]   # will be a DataFrame
Df['lower']  # will be a Series, the distinction is that Series doesn' t have columns

Df.iloc[1]  # will be a series object
Df.iloc[0:2]  # will be a DataFrame
Df.iloc[[0,2]] # will bea DataFrame


Df.iloc[[0,2],[0,2]]   # will be a DataFrame
Df.iloc[1,1]    # will not be a string, the type of each entry itself
Df.iloc[:,0:1]   # will ba a DataFrame, not a Series

Df.iloc[1][['a','b']]   # numeric for row, string for column
Df[['a','b']].iloc[1]  # another way to do that.

Df['new']=[4,5,6]
Df.insert(1,'old',[5,6,7])  # new column will become column 1






df = pd.read_csv('mRNA-ExonIDs.txt',sep='\t',header=None,names=['EnsGID','EnsTID','EnsPID','Exons'])
df = pd.read_csv('GTEx_Analysis_2017-06-05_v8_RSEMv1.3.0_transcript_tpm.gct',sep='\t',skiprows=2,usecols=['transcript_id','gene_id'])

