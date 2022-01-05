import numpy as np
import pandas as pd
import os,sys
import dash
from dash import dcc,html,dash_table
import plotly.graph_objects as go
from dash.dependencies import Input,Output,State
import networkx as nx
from functools import lru_cache
from tqdm import tqdm


# read the data
os.chdir('/Users/ligk2e/Desktop/tmp')
df = pd.read_csv('input.csv',index_col=0)
df['id'] = ['_'.join([s,t,c]) for s,t,c in zip(df['source'],df['timepoint'],df['cell_type'])]

# calculate overlap and build the graph
@lru_cache()
def calculate(id1,id2):
    df1 = df.loc[df['id']==id1,:]
    df2 = df.loc[df['id']==id2,:]
    o = list(set(df1.index).intersection(set(df2.index)))
    lo = len(o)   # length of overlap
    so = ','.join(o)  # stringify the overlap list
    return lo,so

all_ids = df['id'].unique()
tuple_list = []
size_dict = {}
for id1 in tqdm(all_ids):
    size_dict[id1] = df.loc[df['id']==id1,:].shape[0]
    for id2 in all_ids:
        if id1 != id2:
            lo,so = calculate(id1,id2)
            tuple_list.append((id1,id2,lo,so))


pdf = pd.DataFrame.from_records(tuple_list,columns=['source','target','weight','attr'])
G = nx.from_pandas_edgelist(pdf,source='source',target='target',edge_attr=True)
nx.set_node_attributes(G,nx.circular_layout(G),'coord')
nx.set_node_attributes(G,size_dict,'size')
all_nodes = list(G.nodes)


# dash app
app = dash.Dash(__name__)

@app.callback(
    Output('network','figure'),
    Input('select_node','value')
)
def network_graph(selected_nodes):
    sub_G = G.subgraph(nodes=selected_nodes)
    # first draw lines and middle nodes
    edge_x = []
    edge_y = []
    middle_x = []
    middle_y = []
    middle_text = []
    middle_customdata = []
    for edge in sub_G.edges(data=True):
        x0,y0 = sub_G.nodes[edge[0]]['coord']
        x1,y1 = sub_G.nodes[edge[1]]['coord']
        edge_x.append(x0)
        edge_y.append(y0)
        edge_x.append(x1)
        edge_y.append(y1)
        edge_x.append(None)
        edge_y.append(None)
        middle_x.append((x0 + x1)/2)
        middle_y.append((y0 + y1)/2)
        middle_text.append(edge[2]['attr'])
        middle_customdata.append(','.join([edge[0],edge[1]]))
    edge_trace = go.Scatter(x=edge_x,y=edge_y,mode='lines',line={'width':0.1,'color':'black'})
    middle_node_trace = go.Scatter(x=middle_x,y=middle_y,text=middle_text,hoverinfo='text',customdata=middle_customdata,
                                     mode='markers',marker={'size':5,'color':'green','opacity':1})
    # then draw nodes
    node_x = []
    node_y = []
    node_text = []
    for node in sub_G.nodes(data=True):
        x,y = sub_G.nodes[node[0]]['coord']
        node_x.append(x)
        node_y.append(y)
        node_text.append(node[0])
    node_trace = go.Scatter(x=node_x,y=node_y,mode='markers',text=node_text,hoverinfo='text',marker={'color':'red','size':20})
    # now the figure
    fig = go.Figure(data=[edge_trace,middle_node_trace,node_trace],layout=go.Layout(showlegend=False))
    return fig

@app.callback(
    Output('df','data'),
    Input('network','clickData')
)
def render_df(clickData):
    genes = clickData['points'][0]['text'].split(',')
    s,t = clickData['points'][0]['customdata'].split(',')
    sub_df = df.loc[df['id']==s,:].loc[genes,:].reset_index()
    return sub_df.to_dict(orient='records')


app.layout = html.Div([
    # first part, selection
    html.Div([
        dcc.Dropdown(
            id='select_node',
            options=[{'label':node,'value':node} for node in all_nodes],
            value=[all_nodes[0],all_nodes[1],all_nodes[2]],
            multi=True
        )
    ]),
    # upper part, the network graph
    html.Div([
        dcc.Graph(id='network')
    ]),
    # lower part, the table
    html.Div([
        dash_table.DataTable(id='df',
                             columns=[{'name':i,'id':i} for i in df.reset_index().columns])
    ])
])

app.run_server()







