import numpy as np
import pandas as pd
from sklearn.decomposition import PCA

#get_ipython().magic('matplotlib inline')
import matplotlib.pyplot as plt
import seaborn as sns

from scipy.stats import ks_2samp

import copy

import statistics as st
import sklearn

geno = pd.read_pickle("./output/replication/our/batches1-4_merged/batch1234_geno.p")
pheno = pd.read_pickle("./output/replication/our/batches1-4_merged/batch1234_pheno.p")

import sklearn.tree as tr
from sklearn.datasets import make_classification
from sklearn.ensemble import ExtraTreesClassifier
y = pheno["cbmi"]
def label_to_int(x):
    if x == "lean":
        return 0
    else:
        return 1
y = y.apply(label_to_int).values
x = geno.as_matrix()

from scipy.stats import rankdata as rd

# Build a forest and compute the feature importances
runs = []
genes = list(geno.columns)
#ranks1 = {gene:[] for gene in genes}
ranks1 = {gene:0 for gene in genes}
for i in range(50,60):
    forest = ExtraTreesClassifier(n_estimators=10000,
                                  random_state=i)
    
    forest.fit(x, y)

    importances = forest.feature_importances_
    std = np.std([tree.feature_importances_ for tree in forest.estimators_],
                 axis=0)
    indices = np.argsort(importances)[::-1]
     
    #t=rd(importances, method='average')
    #for l in range(0,len(t)):
    #    ranks1[genes[l]].append( t[l])
    
    for l in range(0,len(importances)):
        ranks1[genes[l]]+=importances[l]
        

#for gene in ranks1.keys():
#    ranks1[gene] = st.median(ranks1[gene])


import operator
sorted_genes1 = sorted(ranks1.items(), key=operator.itemgetter(1))[::-1]
signature = pd.DataFrame(sorted_genes1[:100])
signature = signature.set_index(0)

'''
#ranks2 = {gene:[] for gene in genes}
ranks2 = {gene:0 for gene in genes}

for i in range(60,70):
    forest = ExtraTreesClassifier(n_estimators=2500,
                                  random_state=i)
    
    forest.fit(x, y)
    importances = forest.feature_importances_
    std = np.std([tree.feature_importances_ for tree in forest.estimators_],
                 axis=0)
    indices = np.argsort(importances)[::-1]
    
#    t=rd(importances, method='average')
#    for l in range(0,len(t)):
#        ranks2[genes[l]].append(t[l])
    
    for l in range(0,len(importances)):
        ranks2[genes[l]]+=importances[l]

#for gene in ranks2.keys():
#    ranks2[gene] = st.median(ranks2[gene])

    
sorted_genes2 = sorted(ranks2.items(), key=operator.itemgetter(1))[::-1]


s1 = set([sorted_genes1[p][0][:200] for p in range(0,200)])
s2 = set([sorted_genes2[p][0][:200] for p in range(0,200)])
k=list(set.intersection(s1,s2))
'''


# ### Save signature
# We save the signature in pickled and CSV form, as well as a list of all the genes in the final gene-expression matrix. These are known as the **background** genes, and must be used to correctly identify overrepresented gene sets.

# In[22]:
hugo_df = pd.read_table("./data/HUGO_official_list_20160613.txt")
entrez_to_genesymb = dict(hugo_df[["Entrez Gene ID","Approved Symbol"]].astype("unicode").applymap(lambda x:x.split(".")[0]).values)


# In[13]:

signature["GeneSymb"] = signature.index.map(lambda x:entrez_to_genesymb[x] if x in entrez_to_genesymb.keys() else "???")
signature.index.name = "Entrez"
signature.set_index("GeneSymb",append=True,inplace=True)
signature.to_csv("./output/results/randomForest_signature.csv", header = None)
signature.to_pickle("./output/results/randomForest_signature.p")

# In[23]:

#np.savetxt("./output/replication/our/batches1-4_merged/background_genes.txt",v_df.index.get_level_values(0).values,fmt="%s")
