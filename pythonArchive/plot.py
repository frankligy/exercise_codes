#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Apr 17 16:57:31 2020

@author: ligk2e
"""
import numpy as np
import matplotlib.pyplot as plt

# use a scatter plot to illustrate lots of features of a single figure
x = np.linspace(0,100,100)
y = [i+2 for i in x]
palette = ['#BF5842','#AA3B23','#6B5651','#9926A3','#26A348'] 

fig = plt.figure()          # figure object, it is the canvas
plt.scatter(x,y,s=50,c=palette[0],marker = 'o',edgecolor = 'black',label = 'xvalue')
plt.xlabel('this is x-axis')
plt.ylabel('this is y-axis')
plt.title('this is the title')
plt.legend(bbox_to_anchor=(1.04,1),fontsize = 10)
plt.xticks([0.0,20,30,40,50,60],[1,2,3,4,5,6],rotation = 20)
plt.tick_params(axis='both',which='both',bottom=True,labelbottom=True)  # bottom means ticks, label bottom means ticklabels
plt.grid()
plt.text(50,50,'center position')
plt.hlines([75,80],30,50)
plt.savefig('text.pdf',bbox_inches='tight')


# figure object is different from axe object, they might have diverse method, so let's introduce axes
fig = plt.figure()
ax = fig.axes([0,0,1,1]) # this ax anchor at [0,0], width is 100% of canvas, height is 100% of canvas. So it is same as fig
# then we could use ax as handle to do any operations, enjoy!


# subplot
fig,axes = plt.subplots(nrows=3,ncols=4)   # fig is still the whole canvas, axes is a matrix, each entry represents a subplot
axes[2,2].set_title('hi')   # now it is a AxesSubplot object

fig = plt.figure()
ax1 = fig.add_subplot(2,4,1)   # now it is a AxesSubplot object as well





# after the basic of plotting in python, here is the advanced
# figure out what information each kind of plot could convey to the reader?
# catalogue1: barchart, boxchart, violinchart
# catalogue2: histogram to illustrate distribution
# catalogue3: scatter plot
# catalogue4: line plot ......

# for each single plot type, figure out the basic building block and the process of stacking all the features in it

# read the documents, and other resources
# https://python-graph-gallery.com  

# remember the sample plot will help us



# example1: 
# bar chart(hstack), only two elements
labels = ['G1','G2','G3','G4','G5']
men = [45,78,48,78,90]
women = [34,56,89,23,60]

x = np.arange(len(labels))
width = 0.35

fig = plt.figure()
ax = fig.add_axes([0,0,1,1])
ax.bar(x-width/2,men,width,label='men')   # could assign this command to a variable, then this variable could be handled later
ax.bar(x+width/2,women,width,label='women')
ax.set_ylabel('value')
ax.set_xlabel('groups')
ax.set_xticks(x)
ax.set_xticklabels(labels)
ax.legend()
fig.tight_layout()
plt.show()

# example2: 
# bar chart(htack), multiple elements (context please refer to protein_secondary_structure_predictor.py)
import matplotlib.pyplot as plt

fig = plt.figure()

barWidth = 0.9
# in following: 1 means l=5, 2 means l=7, 3 means l=9, 4 means l=11
r1 = [1,5,9,13]
r2 = [2,6,10,14]
r3 = [3,7,11,15]
r4 = [4,8,12,16]
r5 = sorted(r1 + r2 + r3 + r4)

bar1 = [confidenceInterval(item)[0] for item in [accuracy_l5_k3,accuracy_l5_k5,accuracy_l5_k7,accuracy_l5_k10]]
bar2 = [confidenceInterval(item)[0] for item in [accuracy_l7_k3,accuracy_l7_k5,accuracy_l7_k7,accuracy_l7_k10]]
bar3 = [confidenceInterval(item)[0] for item in [accuracy_l9_k3,accuracy_l9_k5,accuracy_l9_k7,accuracy_l9_k10]]
bar4 = [confidenceInterval(item)[0] for item in [accuracy_l11_k3,accuracy_l11_k5,accuracy_l11_k7,accuracy_l11_k10]]

yer1 = np.transpose(np.array([confidenceInterval(item)[2] for item in [accuracy_l5_k3,accuracy_l5_k5,accuracy_l5_k7,accuracy_l5_k10]]))
yer2 = np.transpose(np.array([confidenceInterval(item)[2] for item in [accuracy_l7_k3,accuracy_l7_k5,accuracy_l7_k7,accuracy_l7_k10]]))
yer3 = np.transpose(np.array([confidenceInterval(item)[2] for item in [accuracy_l9_k3,accuracy_l9_k5,accuracy_l9_k7,accuracy_l9_k10]]))
yer4 = np.transpose(np.array([confidenceInterval(item)[2] for item in [accuracy_l11_k3,accuracy_l11_k5,accuracy_l11_k7,accuracy_l11_k10]]))    

plt.bar(r1,bar1,width=barWidth,color=(0.3,0.1,0.4,0.6),yerr = yer1,capsize=4,label='window length=5')
plt.bar(r2,bar2,width=barWidth,color=(0.3,0.33,0.4,0.6),yerr = yer2,capsize=4,label='window length=7')
plt.bar(r3,bar3,width=barWidth,color=(0.3,0.65,0.4,0.6),yerr = yer3,capsize=4,label='window length=9')
plt.bar(r4,bar4,width=barWidth,color=(0.3,0.9,0.4,0.6),yerr = yer4,capsize=4,label='window length=11')

plt.vlines(4.50,0.0,0.60,linestyles='dashed')
plt.vlines(8.50,0.0,0.60,linestyles='dashed')
plt.vlines(12.50,0.0,0.60,linestyles='dashed')
text = ['k=3','k=5','k=7','k=10']
for i in range(4):
    plt.text(x=i*4+2.0,y=0.58,s=text[i],size=12)

plt.legend(bbox_to_anchor=(1.04,1),fontsize=10)
plt.xticks([r for r in r5],['k=3,l=5','k=3,l=7','k=3,l=9','k=3,l=11',
                                       'k=5,l=5','k=5,l=7','k=5,l=9','k=5,l=11',
                                       'k=7,l=5','k=7,l=7','k=7,l=9','k=7,l=11',
                                       'k=10,l=5','k=10,l=7','k=10,l=9','k=10,l=11'],rotation=60)
plt.show()

# svg
import pygal
import cairosvg








