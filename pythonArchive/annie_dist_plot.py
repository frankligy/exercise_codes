import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import matplotlib as mpl

mpl.use('Agg')

count_cite = pd.read_csv('/Users/ligk2e/Downloads/AS_CITE_HSC_ADT_clean.txt',sep='\t',index_col=0)
count_tea = pd.read_csv('/Users/ligk2e/Downloads/h1_TEA_filtered_adt_mtx.txt',sep='\t',index_col=0)

valid_cite = [item for item in count_cite.index if 'isotype_Ctrl' not in item]
valid_tea = [item for item in count_tea.index if 'isotype_Ctrl' not in item]
common = list(set(valid_cite).intersection(set(valid_tea)))

valid_count_cite = count_cite.loc[common,:]  # [108 rows x 9021 columns]
valid_count_tea = count_tea.loc[common,:]   # [108 rows x 6628 columns]

def plot_running_sum(ax,row_cite,row_tea):
    # row_cite and row_tea are both pandas series
    total_cite, total_tea = row_cite.sum(),row_tea.sum()
    cite_hist,_ = np.histogram(row_cite.values,bins=np.logspace(0,4,401))
    tea_hist,_ = np.histogram(row_tea.values,bins=np.logspace(0,4,401))
    data = np.vstack((tea_hist/total_tea,cite_hist/total_cite))
    ax.imshow(data,aspect='auto',interpolation='none',vmin=0.000,vmax=0.001)
    ax.set_yticks([0,1])
    ax.set_yticklabels(['tea','cite'],fontsize=3)
    ax.set_xticks([0,100,200,300,400])
    ax.tick_params(axis='both',length=0.3)
    ax.set_xticklabels(['0','10','100','1k','10k'],fontsize=3)
    ax.set_title(row_cite.name,fontsize=3,color='red',fontweight='bold',pad=0.4)
    ax.axhline(y=0.5,color='w',linewidth=1)
    for direction in ['right','bottom','left','top']:
        ax.spines[direction].set_visible(False)


fig,axes = plt.subplots(nrows=27,ncols=4,figsize=(10,20),gridspec_kw={'wspace':0.5,'hspace':0.5})
axes = axes.flatten()
for i in range(len(axes)):
    ax = axes[i]
    row_cite = valid_count_cite.iloc[i,:]
    row_tea = valid_count_tea.iloc[i,:]
    plot_running_sum(ax,row_cite,row_tea)
plt.savefig('plot.pdf',bbox_inches='tight')
plt.close()






