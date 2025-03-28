#!/gpfs/data/yarmarkovichlab/Frank/immunopeptidome_project/engine/SNV/snv_env/bin/python3.7

import pandas as pd
import numpy as np
import sys,os
from Bio.SeqIO.FastaIO import SimpleFastaParser
from Bio import SeqIO
import subprocess
import matplotlib.pyplot as plt
import matplotlib as mpl
from tqdm import tqdm

mpl.rcParams['pdf.fonttype'] = 42
mpl.rcParams['ps.fonttype'] = 42
mpl.rcParams['font.family'] = 'Arial'


import pycircos
Garc = pycircos.Garc
Gcircle = pycircos.Gcircle

# set up plus tickplot
circle = Gcircle(figsize=(8,8))
with open('example/example_data_chromosome_general.csv') as f:
    f.readline()
    for line in f:
        line = line.rstrip().split(',')
        name = line[0]
        length = int(line[-1])
        arc = Garc(arc_id=name,size=length,interspace=2,raxis_range=(935,985),labelposition=80,label_visible=True)
        circle.add_garc(arc)
circle.set_garcs(-65,245) 
for arc_id in circle.garc_dict:
    circle.tickplot(arc_id, raxis_range=(985,1000), tickinterval=20000000, ticklabels=None) 


# cytoband
color_dict   = {"gneg":"#FFFFFF00", "gpos25":"#EEEEEE", "gpos50":"#BBBBBB", "gpos75":"#777777", "gpos100":"#000000", "gvar":"#FFFFFF00", "stalk":"#C01E27", "acen":"#D82322"}
arcdata_dict = {}
with open('example/example_data_chromosome_cytoband.csv') as f:
    f.readline()
    for line in f:
        line = line.rstrip().split(',')
        name = line[0]
        start = int(line[1]) - 1
        width = int(line[2]) - start
        if name not in arcdata_dict:
            arcdata_dict[name] = {}
            arcdata_dict[name]['positions'] = []
            arcdata_dict[name]['widths'] = []
            arcdata_dict[name]['colors'] = []
        arcdata_dict[name]['positions'].append(start)
        arcdata_dict[name]['widths'].append(width)
        arcdata_dict[name]['colors'].append(color_dict[line[-1]])
for key in arcdata_dict:
    circle.barplot(key,data=[1]*len(arcdata_dict[key]['positions']),positions=arcdata_dict[key]['positions'],
                   width=arcdata_dict[key]['widths'],raxis_range=[935,985],facecolor=arcdata_dict[key]['colors'])

# scatter plot
values_all = []
arcdata_dict = {}
with open('example/example_data_point.csv') as f:
    f.readline()
    for line in f:
        line = line.rstrip().split(',')
        name = line[0]
        start = int(line[1]) - 1
        end = int(line[2])
        mid = (start + end) / 2
        value = float(line[-1])
        values_all.append(value)
        if name not in arcdata_dict:
            arcdata_dict[name] = {}
            arcdata_dict[name]['positions'] = []
            arcdata_dict[name]['values'] = []
        arcdata_dict[name]['positions'].append(mid)
        arcdata_dict[name]['values'].append(value)
vmin,vmax = min(values_all),max(values_all)
for key in arcdata_dict:
    circle.scatterplot(key,data=arcdata_dict[key]['values'],positions=arcdata_dict[key]['positions'],
                       rlim=[vmin-0.05*abs(vmin),vmax+0.05*abs(vmax)],raxis_range=(845,925),facecolor='orange',spine=True)
 
# lineplot
values_all   = [] 
arcdata_dict = {}
with open("example/example_data_point.csv") as f:
    f.readline()
    for line in f:
        line  = line.rstrip().split(",")
        name  = line[0]     
        start = int(line[1])-1
        end   = int(line[2]) 
        mid   = (start+end)/2
        value = float(line[-1]) 
        values_all.append(value) 
        if name not in arcdata_dict:
            arcdata_dict[name] = {}
            arcdata_dict[name]["positions"] = []
            arcdata_dict[name]["values"] = []
        arcdata_dict[name]["positions"].append(mid) 
        arcdata_dict[name]["values"].append(value)
    
vmin, vmax = min(values_all), max(values_all) 
for key in arcdata_dict:
    circle.lineplot(key, data=arcdata_dict[key]["values"], positions=arcdata_dict[key]["positions"], 
                    rlim=[vmin-0.05*abs(vmin), vmax+0.05*abs(vmax)], raxis_range=(755,835), linecolor="royalblue", spine=False)

#bar plot
values_all   = [] 
arcdata_dict = {}
with open("example/example_data_point.csv") as f:
    f.readline()
    for line in f:
        line  = line.rstrip().split(",")
        name  = line[0]     
        start = int(line[1])-1
        end   = int(line[2]) 
        width = end-start 
        if name not in arcdata_dict:
            arcdata_dict[name] = {}
            arcdata_dict[name]["positions"] = []
            arcdata_dict[name]["widths"]    = [] 
            arcdata_dict[name]["values"]    = [] 
        arcdata_dict[name]["positions"].append(start) 
        arcdata_dict[name]["widths"].append(width)
        arcdata_dict[name]["values"].append(float(line[-1]))
        values_all.append(float(line[-1]))

vmin, vmax = min(values_all), max(values_all) 
for key in arcdata_dict:  
    circle.barplot(key, data=arcdata_dict[key]["values"], positions=arcdata_dict[key]["positions"], 
                   width=arcdata_dict[key]["widths"], base_value=0.0, rlim=[vmin-0.05*abs(vmin), vmax+0.05*abs(vmax)],
                   raxis_range=[665,745], facecolor="y", spine=True)

#heatmap
values_all   = [] 
arcdata_dict = {}
with open("example/example_data_point.csv") as f:
    f.readline()
    for line in f:
        line  = line.rstrip().split(",")
        name  = line[0]     
        start = int(line[1])-1
        end   = int(line[2]) 
        width = end-start 
        if name not in arcdata_dict:
            arcdata_dict[name] = {}
            arcdata_dict[name]["positions"] = []
            arcdata_dict[name]["widths"]    = [] 
            arcdata_dict[name]["values"]    = [] 
        arcdata_dict[name]["positions"].append(start) 
        arcdata_dict[name]["widths"].append(width)
        arcdata_dict[name]["values"].append(float(line[-1]))
        values_all.append(float(line[-1]))

vmin, vmax = min(values_all), max(values_all) 
for key in arcdata_dict:
    circle.heatmap(key, data=arcdata_dict[key]["values"], positions=arcdata_dict[key]["positions"], 
                   width=arcdata_dict[key]["widths"], raxis_range=[615,655], vmin=vmin, vmax=vmax, 
                   cmap=plt.cm.viridis)

#linkplot
with open("example/example_data_links.csv") as f:
    f.readline()
    for line in f:
        line  = line.rstrip().split(",")
        name1  = line[0]     
        start1 = int(line[1])-1
        end1   = int(line[2])
        name2  = line[3]     
        start2 = int(line[4])-1
        end2   = int(line[5])
        source = (name1, start1, end1, 615)
        destination = (name2, start2, end2, 615)
        circle.chord_plot(source, destination, facecolor=circle.garc_dict[name1].facecolor)

circle.save('check6')




'''
Go to NCBI assembly to download the species reference assembly, including
genome fasta
annotation gtf
cds fasta
protein fasta

But there are 35 assembly for different strains

color code the cds

only look for ones that mapped to reference assembly
'''



genome_dict = {}
with open('GCA_013267435.1/GCA_013267435.1_ASM1326743v1_genomic.fna','r') as in_handle:
    for title,seq in SimpleFastaParser(in_handle):
        genome_dict[title] = seq

tmp = subprocess.run('grep -A1 \'kraken:taxid|1397\' TCGA-25-1877-01A/test_cseqs_1.fq',shell=True,stdout=subprocess.PIPE,universal_newlines=True).stdout.split('\n')[:-1]
seqs = [item for item in tmp if item != '--' and 'kraken' not in item]

ref = genome_dict['CP053989.1 Niallia circulans strain FDAARGOS_783 chromosome, complete genome']
true_seqs = [item for item in seqs if item in ref]

pos = []
for item in true_seqs:
    i = ref.find(item)
    pos.append(i)

print(len(ref));sys.exit('stop')

fig,ax = plt.subplots()
ax.bar(x=np.arange(len(ref)),y=[1 if i in pos else 0 for i in np.arange(len(ref))])
plt.savefig('check.pdf',bbox_inches='tight')
plt.close()

        