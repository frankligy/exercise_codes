#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sun Jan 26 13:46:01 2020

@author: ligk2e
"""

# -*- coding: utf-8 -*-
"""
Spyder Editor

This is a temporary script file.
"""
import pandas as pd
import os
import urllib3
import io
import Bio
import xmltodict
from Bio.Seq import Seq
from Bio.Alphabet import generic_dna, generic_rna


def open_file(path,delimiter):
    df = pd.read_csv(path,sep=delimiter)
    return df

def UID(df, i):
    uid = df['UID'][i]
    #uid = uid.split(':')[1]
    gene = uid.split(':')[0]
    dict = {}
    gene = gene + ':' + uid.split('|')[1].split(':')[0]
    x = uid.split('|')
    dict[gene] = [x[0].split(':')[2]]
    dict[gene].append(x[1].split(':')[1])
    #print ((uid.split('|')[0],gene + ':' + uid.split('|')[1]))
    #return ((uid.split('|')[0],gene + ':' + uid.split('|')[1]))
    #print dict
    #{'gene:ENSid':[E22-33,E34-56]}
    return dict

def get_coordinates(df, i):
    coordinates = df['Coordinates'][i]
    coordinates = coordinates.split('|')
    coords = list()
    for x in coordinates:
        c = x.split(':')[1]
        coords.append((x.split(':')[0],c.split('-')[0],c.split('-')[1]))
    #print coords
    #coords = coords[0],coords[1]
    #[(chr17,former coordinate,latter coordinate),(chr17,former2 coordinate, latter2 coordinate)]
    return coords

def web_scraping(command,url):
    http = urllib3.PoolManager()
    r = http.request(command,url,preload_content=False)
    r.auto_close = False
    buffer = ''
    for line in io.TextIOWrapper(r):
        buffer += str(line)
  
    doc = xmltodict.parse(buffer)
    frag = doc[u'DASDNA'][u'SEQUENCE'][u'DNA']['#text']
    return frag
    
def extract_dna_sequence(coordinates,i):
    #print coordinates
    coordinates = coordinates[i]
    #print(coordinates)
    if int(coordinates[2]) > int(coordinates[1]):
        # (+) strand
        c1 = int(coordinates[1]) - 15
        c2 = coordinates[1]
        c3 = coordinates[2]
        c4 = int(coordinates[2]) + 15
    else:
        # (-) strand
        c3 = int(coordinates[1])
        c4 = int(coordinates[1]) + 15
        c1 = int(coordinates[2]) - 15
        c2 = int(coordinates[2])
    path1 = 'http://genome.ucsc.edu/cgi-bin/das/hg38/dna?segment=' + coordinates[0] + ':'
    path1 = path1 + str(c1) + ',' + str(c2)
    path2 = 'http://genome.ucsc.edu/cgi-bin/das/hg38/dna?segment=' + coordinates[0] + ':'
    path2 = path2 + str(c3) + ',' + str(c4)
    
    frag1 = web_scraping('Get',path1)
    frag2 = web_scraping('Get',path2)
    #file = urllib2.urlopen('http://genome.ucsc.edu/cgi-bin/das/hg19/dna?segment=chr9:82242348,82242363')
    #file = urllib2.urlopen(path1)
#    http = urllib3.PoolManager()
#    r = http.request('Get',path1)
#    data = file.read()
#    file.close()
#    data = xmltodict.parse(data)
#    frag1 = data[u'DASDNA'][u'SEQUENCE'][u'DNA']['#text']
#
#    file = urllib2.urlopen(path2)
#    data = file.read()
#    file.close()
#    data = xmltodict.parse(data)
#    frag2 = data[u'DASDNA'][u'SEQUENCE'][u'DNA']['#text']
    if int(coordinates[2]) > int(coordinates[1]):
        sequence = frag1 + frag2
    else:
        s = Seq(frag1 + frag2, generic_dna)
        sequence = str(s.reverse_complement())

    #print sequence
    return sequence
    #obtained the 30nt dna seq and make sure all are encoded as coding strand(as opposed to template strand)


def translate(dna):
    dna = Seq(dna,generic_dna)
    prot1 = dna.translate(to_stop=False)
    prot2 = dna[1:].translate(to_stop=False)
    prot3 = dna[2:].translate(to_stop=False)
    #dna = str(dna)
    #mrna = string.replace(dna,'t','u')
    #rna = Seq(mrna,generic_rna)
#    f0 = rna.translate()
#    rna = Seq(mrna[1:],generic_rna)
#    f1 = rna.translate()
#    rna = Seq(mrna[2:],generic_rna)
#    f2 = rna.translate()
    dict = {}
    dict['Frame 1'] = str(prot1)
    dict['Frame 2'] = str(prot2)
    dict['Frame 3'] = str(prot3)
    return dict
    #print f0, f1, f2
    
def translate_junction(path,delimiter):
    df = open_file(path,delimiter)
    uid = []
    coords = []
    for i in range(len(df['UID'])):
        uid.append(UID(df,i))
        coords.append(get_coordinates(df,i))
        #uid: [{},{},{}], coords:[[(),()],[],[]]
    dict_seq_prot = {}
    for j in range(len(uid)):
    #for j in range(1):
        dict = {}
        key = list(uid[j].keys())[0]
        #print key
        for i in range(len(uid[j][key])):
            dna = extract_dna_sequence((coords[j]),i)
            if uid[j][key][i] in dict:
                dict[uid[j][key][i]].append(dna)
            else:
                dict[uid[j][key][i]] = [dna]
            dict[uid[j][key][i]].append(translate(dna))
        #print uid[uid.keys()[j]]
        #dict:{'E34-45':[dna_seq,{three frame}]} all events from one gene
        if key in dict_seq_prot:
            dict_seq_prot[key].append(dict)
        else:
            dict_seq_prot[key] = [dict]
    
    return dict_seq_prot
    #dict_seq_prot:{'gene:ENS':dict}

def make_dict(path):
    d={}
    with open(path,'r') as f1:
        for line in f1:
            items = line.split("\t")
            d[items[0]] = items[1].rstrip("\n")
    #print d[d.keys()[2]]
    return d  

def match(df, df_junctions, translations, d):
    pseq = df["pseq"]
    #print pseq
    isoacc = df["isoacc"]
    unique_acc = df["unique_acc"]
    temp_dict = dict()
    temp_dict2 = dict()
    UID = df_junctions["UID"]
    for key in translations.keys():
        #print(key)
        for m in range(len(translations[key])):
            
            for key2 in translations[key][m].keys():
            #print(key2)
            #isin = False
                for key3 in translations[key][m][key2][1].keys():
                #print(key3)
                    for i in range(len(pseq)):
                    #print(pseq[i])
                    #print translations[key][0][key2][1][key3]
                        if translations[key][m][key2][1][key3] in pseq[i]:
#                        isin = True
#                        print(isin)
                            for j in range(len(UID)):
                                if key in UID[j] and key2 in UID[j]:
                                    #######################split uid, check if it's first or second, store
                                    uid_split = UID[j].split("|")
                                    position = -2
                                    for k in range(len(uid_split)):
                                        if key2 in uid_split[k]:
                                            position = k
                                            #print k
                                    if j in temp_dict.keys():
                                        #temp_dict[j].append(key+":"+key2)
                                        temp_dict[j].append((position,isoacc[i]))
                                    else:
                                        #temp_dict[j] = [key+":"+key2]
                                        temp_dict[j] = [(position,isoacc[i])]
                                    break
                        #####################check sequence protein dbase dict for isoform
            #if isin == False:
                #for key3 in translations[key][0][key2][1].keys():
                    for item in d.keys():
                        if translations[key][m][key2][1][key3] in d[item]:
                            for j in range(len(UID)):
                                if key in UID[j] and key2 in UID[j]:
                                    #######################split uid, check if it's first or second, store
                                    uid_split = UID[j].split("|")
                                    position = -2
                                    for k in range(len(uid_split)):
                                        if key2 in uid_split[k]:
                                            position = k
                                            #print k
                                    if j in temp_dict2.keys():
                                        #temp_dict[j].append(key+":"+key2)
                                        temp_dict2[j].append((position,item))
                                    else:
                                        #temp_dict[j] = [key+":"+key2]
                                        temp_dict2[j] = [(position,item)]
                                    break
                    ###############################################
    #print temp_dict
    lst = [["",""] for i in range(len(UID))]
    #print lst
    for key,value in temp_dict.items():
        #print value[0][1]
        #print value[0][0]
        #print key
        for val in value:
            #print(val)
            if lst[key][int(val[0])] == "":
                lst[key][int(val[0])] = str(val[1])
            else:
                lst[key][int(val[0])] = lst[key][int(val[0])] + ", " + str(val[1])
    lst2 = list(range(len(lst)))
    lst3 = list(range(len(lst)))
    lst4 = list(range(len(lst)))
    lst5 = list(range(len(lst)))
    #print lst
    for i in range(len(lst)):
        #print lst[i]
        lst2[i] = lst[i][0]
        lst3[i] = lst[i][1]
    lst = [["",""] for i in range(len(UID))]
    for key, value in temp_dict2.items():
        for val in value:
            if lst[key][int(val[0])] == "":
                lst[key][int(val[0])] = str(val[1])
            else:
                lst[key][int(val[0])] = lst[key][int(val[0])] + ", " + str(val[1])
    for i in range(len(lst)):
        lst4[i] = lst[i][0]
        lst5[i] = lst[i][1]
    #print lst
    #df = pandas.DataFrame.from_dict(df)
    #df['UID_matched'] = lst
    df_junctions['isoacc1'] = lst2
    df_junctions['isoacc2'] = lst3
    df_junctions['ens_acc1'] = lst4
    df_junctions['ens_acc2'] = lst5
    return df_junctions



if __name__ == '__main__':
    translations = translate_junction('/Users/ligk2e/Desktop/TF_events2.txt','\t')
    df = open_file('/Users/ligk2e/Desktop/CCSB_TFIso_Clones_for_Nathan.txt','\t')
    df_junctions = open_file('/Users/ligk2e/Desktop/TF_events2.txt','\t')
    d = make_dict('/Users/ligk2e/Desktop/SEQUENCE-protein-dbase_exoncomp.txt')
    df_junctions_v2 = match(df, df_junctions, translations, d)
    #df_junctions_v2.to_csv('/Users/ligk2e/Desktop/result1.txt',sep='\t',header=True,index=False)
    
   
    
    

    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    