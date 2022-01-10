import pandas as pd
import numpy as np
from colors import pick_n_colors,infer_to_256,to_rgb,colors_for_set
import os,sys
from matplotlib.patches import Rectangle,Patch,PathPatch
import matplotlib.pyplot as plt
from matplotlib.path import Path

col1 = 'pre_TNC'
col2 = 'post_TNC'

df_raw = pd.read_csv('cellfreq_table.csv',index_col=0)
df = df_raw.loc[:,[col1,col2]]
color_dict = colors_for_set(df.index.tolist())
total = df.sum(axis=0).values
fig,ax = plt.subplots()
ax.set_xlim([-0.2,1.2])
ax.set_ylim([-0.2,1.2])
y1,y2 = 1,1
for row in df.itertuples():
    cluster = row.Index
    color = color_dict[cluster]
    tp1 = eval('row.{}'.format(col1)) / total[0]
    tp2 = eval('row.{}'.format(col2)) / total[1]
    # draw block
    rec1 = Rectangle((0.1,y1-tp1),0.3,tp1,edgecolor=color,facecolor=color,lw=1)
    ax.add_patch(rec1)
    rec2 = Rectangle((0.6,y2-tp2),0.3,tp2,edgecolor=color,facecolor=color,lw=1)
    ax.add_patch(rec2)
    # draw line
    start_up = (0.1+0.3,y1)
    end_up = (0.6,y2)
    control_start_up = (0.5,y1)
    control_end_up = (0.5,y2)
    start_down = (0.1+0.3,y1-tp1)
    end_down = (0.6,y2-tp2)
    control_start_down = (0.5,y1-tp1)
    control_end_down = (0.5,y2-tp2)
    p = Path([start_up,control_start_up,control_end_up,end_up,end_down,control_end_down,control_start_down,start_down,start_up],
             [Path.MOVETO,Path.CURVE4,Path.CURVE4,Path.CURVE4,Path.LINETO,Path.CURVE4,Path.CURVE4,Path.CURVE4,Path.CLOSEPOLY])
    pp = PathPatch(p,fc=color,alpha=0.7)
    ax.add_patch(pp)
    y1 = y1 - tp1
    y2 = y2 - tp2
ax.legend(handles=[Patch(color=v) for v in color_dict.values()],labels=[k for k in color_dict.keys()],loc='upper left',bbox_to_anchor=(1,1),ncol=2,frameon=False)
ax.set_xticks([0.25,0.75])
ax.set_xticklabels([col1,col2])
ax.set_yticks([])
ax.spines['top'].set_visible(False)
ax.spines['right'].set_visible(False)
ax.spines['left'].set_visible(False)
plt.savefig('/Users/ligk2e/Downloads/Sankey_plot/{}_{}.pdf'.format(col1,col2),bbox_inches='tight')
plt.close()
sys.exit('stop')






