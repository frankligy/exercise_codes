import sys
import os
import numpy as np
import pandas as pd
from pyteomics import mzml
import matplotlib.pyplot as plt
import matplotlib as mpl
import multiprocessing as mp
import subprocess
from tqdm import tqdm
import json
import argparse 
from dash import Dash, Input, Output, callback, dash_table, html, dcc, State
import re
from ast import literal_eval

mpl.rcParams['pdf.fonttype'] = 42
mpl.rcParams['ps.fonttype'] = 42
mpl.rcParams['font.family'] = 'Arial'




def hla_formatting(pre,pre_type,post_type):
    if pre_type == 'netMHCpan_output' and post_type == 'netMHCpan_input':  # HLA-A*01:01 to HLA-A01:01
        post = [hla.replace('*','') for hla in pre]
    elif pre_type == 'netMHCpan_input' and post_type == 'netMHCpan_output':  # HLA-A01:01 to HLA-A*01:01
        post = [hla[:5] + '*' + hla[5:] for hla in pre]
    elif pre_type == 'netMHCpan_output' and post_type == 'deepimmuno':  # HLA-A*01:01 to HLA-A*0101
        post = [hla.replace(':','') for hla in pre]
    elif pre_type == 'deepimmuno' and post_type == 'netMHCpan_output': # HLA-A*0101 to HLA-A*01:01
        post = []
        for hla in pre:
            first = hla[:6] # HLA-A*
            second = hla[6:8] # 01, it is so far always 2 digit, no case exceed 2 digit
            third = hla[8:] # 01 or 101, it could be 3 digit
            now = first + second + ':' + third
            post.append(now)
    elif pre_type == 'deepimmuno_nostar' and post_type == 'netMHCpan_input':   # HLA-A0101 to HLA-A01:01
        post = []
        for hla in pre:
            first = hla[:5] # HLA-A
            second = hla[5:7] # 01, it is so far always 2 digit, no case exceed 2 digit
            third = hla[7:] # 01 or 101, it could be 3 digit
            now = first + second + ':' + third
            post.append(now)
    return post


@callback(Output('candidate','data'),Input('cancer_dropdown','value'),Input('submit_button','n_clicks'),State('query_peptide','value'))
def filter_table(cancer,n_clicks,query_peptide):

    final = pd.read_csv('/static/{}_final_enhanced.txt'.format(cancer),sep='\t')
    cond = [False if ('[]' in item) and ('(\'HLA-' not in item) else True for item in final['presented_by_each_sample_hla']]
    final = final.loc[cond,:]
    selected_columns = ['pep','typ','source','highest_score','depmap_median','best_pep','n_psm','presented_by_each_sample_hla','additional_query']
    final = final.loc[:,selected_columns]

    if n_clicks == 0:
        return final.to_dict('records')
    peptides = set(peptide.strip() for peptide in query_peptide.split('\n') if peptide.strip())
    filtered_final = final.loc[final['pep'].isin(peptides)]
    return filtered_final.to_dict('records')

    

@callback(Output('output_header','children'),Output('output_text','children'),Input('candidate','active_cell'),Input('candidate', 'data'),Input('candidate','page_current'),Input('candidate','page_size'))
def click_table(active_cell,data,page_current,page_size):
    if active_cell:
        row = active_cell['row'] + page_current * page_size
        antigen_value = data[row]['pep']
        cell_value = data[row][active_cell['column_id']]
        return 'Selected antigen: {}'.format(antigen_value),'Selected cell: {}'.format(cell_value)
    else:
        return 'Selected antigen: None', 'Selected cell: None'


@callback(Output('psm','src'),Input('cancer_dropdown','value'),Input('candidate','active_cell'),Input('candidate', 'data'),Input('candidate','page_current'),Input('candidate','page_size'))
def draw_psm(cancer,active_cell,data,page_current,page_size):
    if active_cell:
        row = active_cell['row'] + page_current * page_size
        cell_value = data[row]['pep']  # must be pep
        pep = cell_value
        if not os.path.exists(os.path.join(assets_dir,'{}_spectrum_{}.png'.format(cancer,pep))):
            print('missing images')
        return app.get_asset_url('{}_spectrum_{}.png'.format(cancer,pep))
   

@callback(Output('differential_1','src'),Output('differential_2','src'),Output('intensity_1','src'),Output('intensity_2','src'),Input('cancer_dropdown','value'),Input('candidate','active_cell'),Input('candidate', 'data'),Input('candidate','page_current'),Input('candidate','page_size'))
def draw_diffential_and_intensity(cancer,active_cell,data,page_current,page_size):
    if active_cell:
        row = active_cell['row'] + page_current * page_size
        cell_value = data[row]['pep']  # must be pep
        pep = cell_value
        typ = data[row]['typ']
        source = data[row]['source']
        source = modify_source_string(source)

        exist_percentile_plot_path = os.path.join(assets_dir,'{}_{}_{}.png'.format(cancer,pep,'percentile'))
        exist_rank_abundance_plot_path = os.path.join(assets_dir,'{}_{}_rank_abundance.png'.format(cancer,pep))

        exist_percentile_plot_file_name = '{}_{}_{}.png'.format(cancer,pep,'percentile')
        exist_rank_abundance_plot_file_name = '{}_{}_rank_abundance.png'.format(cancer,pep)

        

        if typ == 'self_gene':
            ensg,enst,symbol = source.split('|')   
            exist_plot_gene_file_name =  '{}_{}_{}_expr_{}.png'.format(cancer,ensg,symbol,'boxplot+boxplot')
            if (not os.path.exists(os.path.join(assets_dir,exist_plot_gene_file_name))) or (not os.path.exists(exist_percentile_plot_path)) or (not os.path.exists(exist_rank_abundance_plot_path)):
                print('missing images')
            return app.get_asset_url(exist_plot_gene_file_name),None,app.get_asset_url(exist_percentile_plot_file_name),app.get_asset_url(exist_rank_abundance_plot_file_name)

        elif typ == 'splicing':
            coords = source.split('|')[0]
            exist_plot_splicing_file_name = '{}_{}_splicing.png'.format(cancer,coords)
            if (not os.path.exists(os.path.join(assets_dir,exist_plot_splicing_file_name))) or (not os.path.exists(exist_percentile_plot_path)) or (not os.path.exists(exist_rank_abundance_plot_path)):
                print('missing images')
            return app.get_asset_url(exist_plot_splicing_file_name),None,app.get_asset_url(exist_percentile_plot_file_name),app.get_asset_url(exist_rank_abundance_plot_file_name)
        elif typ == 'TE_chimeric_transcript':
            coords = source.split('|')[0]
            erv = source.split('|')[7].split(',')[1]
            exist_plot_splicing_file_name = '{}_{}_splicing.png'.format(cancer,coords)
            exist_plot_erv_file_name = '{}_expr.png'.format(erv)
            if (not os.path.exists(os.path.join(assets_dir,exist_plot_splicing_file_name))) or (not os.path.exists(os.path.join(assets_dir,exist_plot_erv_file_name))) or (not os.path.exists(exist_percentile_plot_path)) or (not os.path.exists(exist_rank_abundance_plot_path)):
                print('missing images')
            return app.get_asset_url(exist_plot_splicing_file_name),app.get_asset_url(exist_plot_erv_file_name),app.get_asset_url(exist_percentile_plot_file_name),app.get_asset_url(exist_rank_abundance_plot_file_name)

        elif typ == 'ERV':
            erv = source.split('|')[0]  
            exist_plot_erv_file_name = '{}_{}_expr.png'.format(cancer,erv)
            if (not os.path.exists(os.path.join(assets_dir,exist_plot_erv_file_name))) or (not os.path.exists(exist_percentile_plot_path)) or (not os.path.exists(exist_rank_abundance_plot_path)):
                print('missing images')
            return  app.get_asset_url(exist_plot_erv_file_name),None,app.get_asset_url(exist_percentile_plot_file_name),app.get_asset_url(exist_rank_abundance_plot_file_name)
        
        elif typ == 'fusion' or typ == 'variant' or typ == 'pathogen' or typ == 'nuORF' or typ == 'intron_retention':  # no diff plot
            if not os.path.exists(exist_percentile_plot_path) or (not os.path.exists(exist_rank_abundance_plot_path)):
                print('missing images')
            return  None,None,app.get_asset_url(exist_percentile_plot_file_name),app.get_asset_url(exist_rank_abundance_plot_file_name)

def modify_source_string(source):
    pat = re.compile(r'ENSG(\d+)\|ENST(\d+)\|')
    if ';' in source:
        sources = source.split(';')
        for item in sources:
            match = re.search(pat,item)
            if match:
                source = item
                break
    if ';' in source:  # further get rid of nc and nuorf
        sources = source.split(';')
        tmp = []
        for item in sources:
            if ('nuORF' not in item) and ('nc|' not in item):
                tmp.append(item)
        source = ';'.join(tmp)

        if ';' in source:
            source = 'not_unique'


    return source


@callback(Output('hla_table','data'),Input('cancer_dropdown','value'),Input('candidate','active_cell'),Input('candidate', 'data'),Input('candidate','page_current'),Input('candidate','page_size'))   
def display_hla_table(cancer,active_cell,data,page_current,page_size):
    meta = pd.read_csv('/static/{}_metadata.txt'.format(cancer),sep='\t')
    if active_cell:
        row = active_cell['row'] + page_current * page_size
        lists = literal_eval(data[row]['additional_query'])
        hlas = literal_eval(data[row]['presented_by_each_sample_hla'])
        for k,vs in hlas.items():
            if len(vs) == 0:
                continue
            else:
                for v in vs:
                    if v[0] is not None:
                        lists.append(v)
        lists = list(set(lists))
        df = pd.DataFrame.from_records(lists,columns=['hla','rank_pert','nM','id'])
        hla = [item.replace('*','') for item in df['hla'].tolist()]
        df['freq'] = [hla2freq.get(item,None) for item in hla]
        # adding recurrency
        hla2info = {}
        considered = {}  # {hla:[s1,s2]}
        for k,vs in hlas.items():
            for v in vs:
                if len(v) == 4:
                    hla,pert,nm,id_ = v
                    if hla is not None:
                        considered.setdefault(hla,[]).append(k)
        lists = literal_eval(data[row]['additional_query'])
        for v in lists:
            hla,pert,nm,id_ = v
            if hla not in considered.keys():
                considered[hla] = []
        for k,vs in considered.items():
            hla = k.split('HLA-')[1]
            sub_meta = meta.loc[meta['HLA'].notna(),:]
            samples = sub_meta.loc[sub_meta['HLA'].str.contains(hla),:]['biology'].unique().tolist()
            try:
                frac = len(vs)/len(samples)
            except ZeroDivisionError:
                frac = 0
            hla2info[k] = (len(vs),len(samples),round(frac,2))
        col1,col2,col3 = [],[],[]
        for hla in df['hla']:
            info = hla2info.get(hla,(0,0,0))
            col1.append(info[0])
            col2.append(info[1])
            col3.append(info[2])
        df['n_detected'] = col1
        df['n_total'] = col2
        df['frac'] = col3
        
        df.sort_values(by='frac',ascending=False,inplace=True)
        
        return df.to_dict('records')



if __name__ == '__main__':

    us_hla = pd.read_csv('/static/US_HLA_frequency.csv',sep=',',index_col=0)
    us_hla.index = hla_formatting(us_hla.index.to_list(),'deepimmuno_nostar','netMHCpan_input')
    hla2freq = us_hla['Percent US population'].to_dict()
    assets_dir = '/assets'

    # start to build app
    app = Dash(__name__,assets_folder='/assets')
    app.layout = html.Div([
        # Main container with two columns
        html.Div([
            # left column, text area, 10%
            html.Div([

                html.H3('Select Cancer'),
                html.H5('Choose one cancer abbreviation',style={'color':'red'}),
                dcc.Dropdown(
                    id = 'cancer_dropdown',
                    options = [
                            'BRCA',
                            'KIRC',
                            'COAD',
                            'STAD',
                            'MESO',
                            'LIHC',
                            'ESCA',
                            'CESC',
                            'BLCA',
                            'RT',
                            'AML',
                            'DLBC',
                            'GBM',
                            'NBL',
                            'PAAD',
                            'HNSC',
                            'OV',
                            'LUSC',
                            'LUAD',
                            'CHOL',
                            'SKCM'
                    ],
                    value = 'NBL',
                    style={'width': '100%'}
                ),


                html.H3('Query Peptide'),
                html.H5('Type in peptides of interest to filter the whole table, one peptide per row please',style={'color':'red'}),
                dcc.Textarea(
                    id='query_peptide',
                    style={'width':'100%','height':'50%','resize':'none'}
                ),
                html.Button('Submit', id='submit_button', n_clicks=0, style={'margin-top': '10px', 'width': '100%'})
            ], style={'width':'10%'}),

            # right column, existing content, 90%
            html.Div([
                html.Div(html.H1('ImmunoVerse Webportal'),style={'text-align':'center'}),
                html.Div(dash_table.DataTable(id='candidate',
                                                page_size=20, page_current=0, page_action='native',
                                                style_data={
                                                    'width': '50px', 'minWidth': '50px', 'maxWidth': '50px',
                                                    'overflow': 'hidden',
                                                    'textOverflow': 'ellipsis',
                                                })),
                html.Div([html.H2(id='output_header'),html.P(id='output_text')]),
                html.Div(html.H2('PSM plot')),
                html.Img(id='psm',width='50%',height='50%',style={'border-style':'dashed'}),

                html.Div([
                        html.H2('Differential plots'),
                        html.Img(id='differential_1',width='45%',height='50%',style={'border-style':'dashed','float':'left'}),
                        html.Img(id='differential_2',width='45%',height='50%',style={'border-style':'dashed','float':'right'})
                    ],style={'overflow':'hidden'}),

                html.Div([
                        html.H2('Intensity plots'),
                        html.Img(id='intensity_1',width='45%',height='50%',style={'border-style':'dashed','float':'left'}),
                        html.Img(id='intensity_2',width='45%',height='50%',style={'border-style':'dashed','float':'right'}),
                    ],style={'overflow':'hidden'}),

                html.Div(html.H2('HLA table')),
                html.Div(dash_table.DataTable(id='hla_table',
                                            page_size=20,page_current=0,page_action='native',   
                                            style_data={
                                                    'width': '50px', 'minWidth': '50px', 'maxWidth': '50px',
                                                    'overflow': 'hidden',
                                                    'textOverflow': 'ellipsis'}))
                
            ],style={'width': '90%', 'padding': '10px', 'box-sizing': 'border-box'})

        ],style={'display': 'flex', 'width': '100vw', 'height': '100vh'})
    ])

    host = subprocess.run(['hostname'],stdout=subprocess.PIPE,universal_newlines=True).stdout.split('\n')[0]
    port = 3838
    app.run(host=host,port=port)