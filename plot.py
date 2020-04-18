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


# bar chart(hstack)
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


# svg
import pygal
import cairosvg








