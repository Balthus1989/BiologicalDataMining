#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Sep 29 11:51:33 2017

@author: Gianmarco Piccinno
"""

import GEOparse
import pandas as pd
import numpy as np
import pylab as pl
import seaborn as sns
import sklearn
from sklearn import *
import math
import subprocess
import os
import re

##################################################################
#               Use only for downloading data!!!
##################################################################

list_data = ["GSE2508", "GSE26637", "GSE27949", "GSE48964"]

gse1 = GEOparse.get_GEO("GSE2508", 
                        destdir="DataMining")
gse1 = GEOparse.get_GEO("GSE26637", 
                        destdir="DataMining")
gse1 = GEOparse.get_GEO("GSE27949", 
                        destdir="DataMining")
gse1 = GEOparse.get_GEO("GSE48964", 
                        destdir="DataMining")

##################################################################
#               Use if data are already available!!!
##################################################################

gse1 = GEOparse.get_GEO(filepath="DataMining/GSE2508_family.soft.gz")
gse2 = GEOparse.get_GEO(filepath="DataMining/GSE26637_family.soft.gz")
gse3 = GEOparse.get_GEO(filepath="DataMining/GSE27949_family.soft.gz")
gse4 = GEOparse.get_GEO(filepath="DataMining/GSE48964_family.soft.gz")

##################################################################
#                           Batch 1!!!
##################################################################

plats_1 = list(gse1.gpls.keys())
plats_1

samples = gse1.phenotype_data[["platform_id", "title"]]; samples
samples.reset_index(level=0, inplace=True); samples
sample = samples.groupby(["platform_id"]); sample.groups

d = {}                        
for l in plats_1:
    print("\nPlatform: "+str(l)+"\n", sample.get_group(l))
    print("\nPlatform: "+str(l)+"\n", sample.get_group(l)['title'])
    ls = "".join(list(sample.get_group(l)['title']))
    lf = re.findall("Lean F", ls)
    of = re.findall("Obese F", ls)
    lm = re.findall("Lean M", ls)
    om = re.findall("Obese M", ls)
    print("LF: ", len(lf), "\nOF: ", len(of), "\nLM: ", len(lm), "\nOM: ", len(om))
    d[l] = {"LF": len(lf), "OF": len(of), "LM": len(lm), "OM": len(om)}
#print(d)
df = pd.DataFrame(d); print(df)
df.sum(axis=1)
df.sum(axis=0)

#preparation of GPL8300/GPL91 platforms ann.
t = sample.get_group("GPL8300")[["index", "title"]]
t = t.set_index('index'); t_l = list(t.index)

with open('GPL8300.txt', 'r+') as handle:
    gse1.gpls['GPL8300'].to_soft(handle)
       

with open('GPL91.txt', 'w') as handle:
    gse1.gpls['GPL91'].to_soft(handle)


#Bash
#sed 's/[#!\^].*$//' < GPL8300.txt > GPL8300w.txt
#sed '/^$/d' GPL8300w.txt > GPL8300wy.txt 

#sed 's/[#!\^].*$//' < GPL91.txt > GPL91w.txt
#sed '/^$/d' GPL91w.txt > GPL91wy.txt

#Bash
#sed 's/[#!\^].*$//' < GSE2508-GPL8300_series_matrix.txt > GSE2508-GPL8300_series_matrix1.txt
#sed '/^$/d' GSE2508-GPL8300_series_matrix1.txt > GSE2508-GPL8300_series_matrix1y.txt 

#sed 's/[#!\^].*$//' < GSE2508-GPL91_series_matrix.txt > GSE2508-GPL91_series_matrix1.txt
#sed '/^$/d' GSE2508-GPL91_series_matrix1.txt > GSE2508-GPL91_series_matrix1y.txt

ann_8300 = pd.read_csv("GPL8300wy.txt", sep="\t", header=0, dtype=str)
ann_8300 = ann_8300[["ID", "Gene Symbol", "ENTREZ_GENE_ID"]]
ann_8300.head()

ann_91 = pd.read_csv("GPL91wy.txt", sep="\t", header=0, dtype=str)
ann_91 = ann_91[["ID", "Gene Symbol", "ENTREZ_GENE_ID"]]
ann_91.head()


ser_8300 = pd.read_csv("GSE2508-GPL8300_series_matrix1y.txt", 
                       sep="\t", header=0, dtype=str); ser_8300.head()

cols_8300 = list(ser_8300.columns)[1:]; cols_8300
d_types_8300 = {}
for el in cols_8300:
    d_types_8300[el] = np.float64

d_types_8300

ser_8300 = pd.read_csv("GSE2508-GPL8300_series_matrix1y.txt", 
                       sep="\t", header=0, dtype=d_types_8300); ser_8300.head()
d_8300 = samples.loc[[el for el in cols_8300]]["title"].to_dict(); d_8300
ser_8300 = ser_8300.rename(index=str, columns=d_8300); ser_8300
                          
ann_ser8300 = ann_8300.merge(ser_8300.rename(columns={'ID_REF':'ID'}),how='inner'); ann_ser8300.head()

#

ser_91 = pd.read_csv("GSE2508-GPL91_series_matrix1y.txt", 
                       sep="\t", header=0, dtype=str); ser_8300.head()

cols_91 = list(ser_91.columns)[1:]; cols_91
d_types_91 = {}
for el in cols_91:
    d_types_91[el] = np.float64

d_types_91

ser_91 = pd.read_csv("GSE2508-GPL91_series_matrix1y.txt", 
                       sep="\t", header=0, dtype=d_types_91); ser_8300.head()

d_91 = samples.loc[[el for el in cols_91]]["title"].to_dict(); d_91
ser_91 = ser_91.rename(index=str, columns=d_91); ser_91

ann_ser91 = ann_91.merge(ser_91.rename(columns={'ID_REF':'ID'}),how='inner'); ann_ser91.head()

##################################################################
#                           Batch 2!!!
##################################################################

plats_2 = list(gse2.gpls.keys())[0]; plats_2

samples = gse2.phenotype_data[["title"]]; samples

                            
with open('GPL570.txt', 'w') as handle:
    gse2.gpls['GPL570'].to_soft(handle)
    
#Bash
#sed 's/[#!\^].*$//' < GPL570.txt > GPL570w.txt
#sed '/^$/d' GPL570w.txt > GPL570wy.txt 

#Bash
#sed 's/[#!\^].*$//' < GSE26637_series_matrix.txt > GSE26637_series_matrix1.txt
#sed '/^$/d' GSE26637_series_matrix1.txt > GSE26637_series_matrix1y.txt

ann_570 = pd.read_csv("GPL570wy.txt", sep="\t", header=0, dtype=str)
ann_570 = ann_570[["ID", "Gene Symbol", "ENTREZ_GENE_ID"]]; ann_570.head()

ser_570 = pd.read_csv("GSE26637_series_matrix1y.txt", 
                       sep="\t", header=0, dtype=str); ser_570.head()

cols_570 = list(ser_570.columns)[1:]; cols_570
d_types_570 = {}
for el in cols_570:
    d_types_570[el] = np.float64

d_types_570

ser_570 = pd.read_csv("GSE26637_series_matrix1y.txt", 
                       sep="\t", header=0, dtype=d_types_570); ser_570.head()

d_570 = samples.loc[[el for el in cols_570]]["title"].to_dict(); d_570
ser_570 = ser_570.rename(index=str, columns=d_570); ser_570
                          
ann_ser570 = ann_570.merge(ser_570.rename(columns={'ID_REF':'ID'}),how='inner'); ann_ser570.head()


##################################################################
#                           Batch 3!!!
##################################################################

plats_3 = list(gse3.gpls.keys())[0]; plats_3

samples = gse3.phenotype_data[["title"]]; samples
                             
with open('GPL570v2.txt', 'w') as handle:
    gse3.gpls['GPL570'].to_soft(handle)

#Bash
#sed 's/[#!\^].*$//' < GPL570v2.txt > GPL570v2w.txt
#sed '/^$/d' GPL570v2w.txt > GPL570v2wy.txt 

#Bash
#sed 's/[#!\^].*$//' < GSE27949_series_matrix.txt > GSE27949_series_matrix1.txt
#sed '/^$/d' GSE27949_series_matrix1.txt > GSE27949_series_matrix1y.txt

ann_570v2 = pd.read_csv("GPL570v2wy.txt", sep="\t", header=0, dtype=str)
ann_570v2 = ann_570v2[["ID", "Gene Symbol", "ENTREZ_GENE_ID"]]; ann_570v2.head()
ser_570v2 = pd.read_csv("GSE27949_series_matrix1y.txt", 
                       sep="\t", header=0, dtype=str); ser_570v2.head()

cols_570v2 = list(ser_570v2.columns)[1:]; cols_570v2
                 
d_types_570v2 = {}
for el in cols_570v2:
    d_types_570v2[el] = np.float64

d_types_570v2

ser_570v2 = pd.read_csv("GSE27949_series_matrix1y.txt", 
                       sep="\t", header=0, dtype=d_types_570v2); ser_570v2.head()

d_570v2 = samples.loc[[el for el in cols_570v2]]["title"].to_dict(); d_570v2

ser_570v2 = ser_570v2.rename(index=str, columns=d_570v2); ser_570v2.head()

ann_ser570v2 = ann_570v2.merge(ser_570v2.rename(columns={'ID_REF':'ID'}),how='inner'); ann_ser570v2.head()

##################################################################
#                           Batch 4!!!
##################################################################

plats_4 = list(gse4.gpls.keys())[0]; plats_4

samples = gse4.phenotype_data[["title"]]; samples

with open('GPL6244.txt', 'w') as handle:
    gse4.gpls['GPL6244'].to_soft(handle)

#Bash
#sed 's/[#!\^].*$//' < GPL6244.txt > GPL6244w.txt
#sed '/^$/d' GPL6244w.txt > GPL6244wy.txt 

#Bash
#sed 's/[#!\^].*$//' < GSE48964_series_matrix.txt > GSE48964_series_matrix1.txt
#sed '/^$/d' GSE48964_series_matrix1.txt > GSE48964_series_matrix1y.txt

ann_GPL6244 = pd.read_csv("GPL6244wy.txt", sep="\t", header=0, dtype=str).dropna()
ann_GPL6244.head()
ann_GPL6244.gene_assignment = ann_GPL6244.gene_assignment.apply(func=lambda x: x.split(sep='//')[0:2])
ann_GPL6244 = ann_GPL6244[["ID", "gene_assignment"]]; ann_GPL6244.head()
ann_GPL6244["gene_symbol"] = ann_GPL6244.gene_assignment.apply(lambda s: s[-1].strip()); ann_GPL6244.head()
ann_6244 = ann_GPL6244[["ID", "gene_symbol"]]; ann_6244.head()  #annotation_table

ser_6244 = pd.read_csv("GSE48964_series_matrix1y.txt", 
                       sep="\t", header=0, dtype=str); ser_6244.head()

cols_6244 = list(ser_6244.columns); cols_6244
                 
d_types_6244 = {}
for el in cols_6244:
    if el == "ID_REF":
        d_types_6244[el] = str
    else:
        d_types_6244[el] = np.float64

d_types_6244

ser_6244 = pd.read_csv("GSE48964_series_matrix1y.txt", 
                       sep="\t", header=0, dtype=d_types_6244); ser_6244.head()

d_6244 = samples.loc[[el for el in cols_6244[1:]]]["title"].to_dict(); d_6244

ser_6244 = ser_6244.rename(index=str, columns=d_6244); ser_6244.head()

ann_ser6244 = ann_6244.merge(ser_6244.rename(columns={'ID_REF':'ID'}),how='inner')
ann_ser6244.head()

ann4 = set(list(ann_ser6244['ID']))
len(list(ann_ser6244['ID']))
len(ann4)

def rem_dups(annotations):
    grouped = annotations.groupby('ID').filter(lambda x: len(set(list(x['gene_symbol'])))==1)
    g_ann = grouped.drop_duplicates()
    g_ann = g_ann[g_ann.gene_symbol.str.contains('[A-Z1-9]')]
    return g_ann

g_ann4 = rem_dups(ann_ser6244)
g_ann4.shape
g_ann4.head()
g_ann4 = g_ann4.groupby('gene_symbol').mean()
g_ann4.shape
g_ann4.head()

with open('symbols.txt', 'w') as handle:
    g_ann4.to_csv(handle, sep='\n', columns=[], header=False, index=True)
    
#upload the file symbols.txt in DAVID and gain the EntrezIDs

#       DATA TRANSFORMATIONS

g_ann4 = g_ann4.apply(lambda x: np.log2(x))
g_ann4.head()
x = g_ann4.values
min_max_scaler = preprocessing.MinMaxScaler()
x_scaled = min_max_scaler.fit_transform(x)
df = pd.DataFrame(x_scaled)
