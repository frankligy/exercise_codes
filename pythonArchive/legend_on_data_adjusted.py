
import matplotlib.pyplot as plt
from adjustText import adjust_text
import pandas as pd
from sctriangulate.colors import *

umap = pd.read_csv('/Users/ligk2e/Downloads/FW__CITE-HIVE_Integration_UMAP_Coordinates/UMAP_Coords_CITE_HIVE.txt',sep='\t',index_col=0)
group = pd.read_csv('/Users/ligk2e/Downloads/FW__CITE-HIVE_Integration_UMAP_Coordinates/groups.r7-10x-hive.txt',sep='\t',index_col=0,header=None)
combined = umap.join(group,how='inner')
c = colors_for_set(list(combined[2].unique()))
combined['color'] = combined[2].map(c).values
cluster_to_centroid = []
for cluster,sub_df in combined.groupby(by=2):
    x_mean = sub_df['UMAP_1'].values.mean()
    y_mean = sub_df['UMAP_2'].values.mean()
    cluster_to_centroid.append((x_mean,y_mean,cluster))
fig,ax = plt.subplots()
ax.scatter(combined['UMAP_1'],combined['UMAP_2'],s=0.5,c=combined['color'])
texts = [ax.text(x=x,y=y,s=s,fontsize=4) for x,y,s in cluster_to_centroid]
adjust_text(texts,arrowprops=dict(arrowstyle='->', color='red'))
plt.savefig('/Users/ligk2e/Desktop/test.pdf',bbox_inches='tight')
plt.close()

import plotly.graph_objects as go
node_traces = []
for cluster,sub_df in combined.groupby(by=2):
    node_trace = go.Scatter(name='nodes',x=sub_df['UMAP_1'].tolist(),y=sub_df['UMAP_2'].tolist(),mode='markers',
                            hoverinfo='text',text=sub_df[2].tolist(),marker={'color':sub_df['color'].iloc[0],'size':1})
    node_traces.append(node_trace)
fig_layout = go.Layout(showlegend=False)
fig = go.Figure(data=node_traces,layout=fig_layout)
fig.write_html('/Users/ligk2e/Desktop/test.html',include_plotlyjs='cdn')
