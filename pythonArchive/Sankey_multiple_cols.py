import pandas as pd
import numpy as np
from colors import pick_n_colors,infer_to_256,to_rgb,colors_for_set
import os,sys
from matplotlib.patches import Rectangle,Patch,PathPatch
import matplotlib.pyplot as plt
from matplotlib.path import Path
import matplotlib as mpl

mpl.rcParams['pdf.fonttype'] = 42
mpl.rcParams['ps.fonttype'] = 42
mpl.rcParams['font.family'] = 'Arial'



col_to_plot = ['preTNC','postTNC','relapseTNC','VC2']
df_raw = pd.read_csv('sankey-4col.txt',index_col=0,sep='\t')
alpha = 0.6

df = df_raw.loc[:,col_to_plot]
color_dict = colors_for_set(df.index.tolist())
total = df.sum(axis=0).values
fig,ax = plt.subplots()
ax.set_xlim([0,1])
ax.set_ylim([-0.2,1.2])
ax.spines['top'].set_visible(False)
ax.spines['right'].set_visible(False)
ax.spines['left'].set_visible(False)
margin = 0.1
n = len(col_to_plot)
ys = [1] * n
dist = (1 - margin * 2)/ (n + n -1)
for row in df.itertuples():
    cluster = row.Index
    color = color_dict[cluster]

    for i,col in enumerate(col_to_plot):
        tp = eval('row.{}'.format(col)) / total[i]
        # draw block
        xy_x = 0.1 + i * 2 * dist
        xy_y = ys[i] - tp
        width = dist
        height = tp
        rec = Rectangle((xy_x,xy_y),width,height,edgecolor=color,facecolor=color,lw=1)
        ax.add_patch(rec)
        # draw line
        if i < n-1:
            start_up = (xy_x+dist,ys[i])
            end_up = (xy_x+dist+dist,ys[i+1])
            control_start_up = (xy_x+dist+dist/2,ys[i])
            control_end_up = (xy_x+dist+dist/2,ys[i+1])
            start_down = (xy_x+dist,ys[i]-tp)
            end_down = (xy_x+dist+dist,ys[i+1]-(eval('row.{}'.format(col_to_plot[i+1])) / total[i+1]))
            control_start_down = (xy_x+dist+dist/2,ys[i]-tp)
            control_end_down = (xy_x+dist+dist/2,ys[i+1]-(eval('row.{}'.format(col_to_plot[i+1])) / total[i+1]))
            p = Path([start_up,control_start_up,control_end_up,end_up,end_down,control_end_down,control_start_down,start_down,start_up],
                     [Path.MOVETO,Path.CURVE4,Path.CURVE4,Path.CURVE4,Path.LINETO,Path.CURVE4,Path.CURVE4,Path.CURVE4,Path.CLOSEPOLY])  
            pp = PathPatch(p,fc=color,alpha=alpha)      
            ax.add_patch(pp)
        ys[i] = ys[i] - tp

ax.legend(handles=[Patch(color=v) for v in color_dict.values()],labels=[k for k in color_dict.keys()],loc='upper left',bbox_to_anchor=(1,1),ncol=2,frameon=False)
ax.set_xticks([0.1 + dist/2 + i * 2 * dist for i in range(n)])
ax.set_xticklabels(col_to_plot)
ax.set_yticks([])
plt.savefig('sankey_like_plot_multiple_cols_{}.pdf'.format('_'.join(col_to_plot)),bbox_inches='tight')
plt.close()






