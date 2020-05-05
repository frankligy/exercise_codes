#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon May  4 16:13:18 2020

@author: ligk2e
"""

import seaborn as sns
import matplotlib.pyplot as plt
import pandas as pd
import numpy as np

df = pd.read_csv('/Users/ligk2e/Downloads/base vs. TA_edited.csv',sep=',')
import random
seq=['1','2']
df['group'] = [ random.choice(seq) for i in range(df.shape[0]) ]

dfgp = df.groupby('group')   # DataFrameGroupBy object
dfgp.first()
group1 = dfgp.get_group('1')['red (base)'].tolist()
group2 = dfgp.get_group('2')['red (base)'].tolist()

group3 = dfgp.get_group('1')['CCT (base)'].tolist()
group4 = dfgp.get_group('2')['CCT (base)'].tolist()

fig = plt.figure()

left, width = 0.1, 0.65
bottom, height = 0.1, 0.65
spacing = 0.005


rect_scatter = [left, bottom, width, height]
rect_histx = [left, bottom + height + spacing, width, 0.2]
rect_histy = [left + width + spacing, bottom, 0.2, height]

ax_scatter = plt.axes(rect_scatter)
ax_scatter.tick_params(direction='in', top=True, right=True)
ax_histx = plt.axes(rect_histx)
ax_histx.tick_params(direction='in', labelbottom=False)
ax_histy = plt.axes(rect_histy)
ax_histy.tick_params(direction='in', labelleft=False)

binwidth=0.50
binsx = np.arange(0,30,binwidth)
binsy = np.arange(5,40,binwidth)

sns.scatterplot(x = "red (base)", y = "CCT (base)", data = df, hue="group",ax=ax_scatter,legend=False,palette=['#0B93F8','#ED8ABF'])
sns.distplot(group1,color='#0B93F8',ax=ax_histx,bins=binsx)
sns.distplot(group2,color='#ED8ABF',ax=ax_histx,bins=binsx)
             
sns.distplot(group3,color='#0B93F8',ax=ax_histy,vertical=True,bins=binsy)
sns.distplot(group4,color='#ED8ABF',ax=ax_histy,vertical=True,bins=binsy)
             
plt.savefig('jointplot.pdf',bbox_inches='tight')