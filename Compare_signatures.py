# -*- coding: utf-8 -*-
"""
Created on Tue Nov 28 12:23:10 2017

@author: damiano
"""
import numpy as np
import pandas as pd

import seaborn as sns


import copy

import statistics as st

sig_rf = pd.read_pickle("./output/results/randomForest_signature2.p")
sig_rf = sig_rf.index.values
sig_gs_rf = []
for gene in sig_rf:
    sig_gs_rf.append(gene[1]) 
sig_gs_rf = sig_gs_rf[:50]
signature = pd.read_pickle("./output/results/signature.p")
signature = signature.index.values
intersection = set.intersection(set (sig_gs_rf), set(signature))
print (len(intersection))

