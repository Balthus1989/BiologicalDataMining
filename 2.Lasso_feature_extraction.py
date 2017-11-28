# -*- coding: utf-8 -*-
"""
Created on Tue Nov 28 14:52:01 2017

@author: damiano
"""
import numpy as np
import pandas as pd

#get_ipython().magic('matplotlib inline')
from sklearn.feature_selection import SelectFromModel
from sklearn.linear_model import LassoCV

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


#Lasso feature selection
clf = LassoCV()

# Set a minimum threshold of 0.25
sfm = SelectFromModel(clf, threshold=0.0000025)
sfm.fit(x, y)
relevant_features_indices = sfm.get_support(indices=True)
print (len(q))

r=sfm.transform(x)
n_features = sfm.transform(x).shape[1]

genes = list(geno.columns)
signature = []
for index in relevant_features_indices:
    signature.append(genes[index])
signature = pd.DataFrame(signature)

# ### Save signature
# We save the signature in pickled and CSV form, as well as a list of all the genes in the final gene-expression matrix. These are known as the **background** genes, and must be used to correctly identify overrepresented gene sets.

# In[22]:
hugo_df = pd.read_table("./data/HUGO_official_list_20160613.txt")
entrez_to_genesymb = dict(hugo_df[["Entrez Gene ID","Approved Symbol"]].astype("unicode").applymap(lambda x:x.split(".")[0]).values)


# In[13]:

signature["GeneSymb"] = signature[0].map(lambda x:entrez_to_genesymb[x] if x in entrez_to_genesymb.keys() else "???")
signature.index.name = "Entrez"
signature.set_index(0,append=False,inplace=True)
signature.to_csv("./output/results/lasso_signature.csv", header = None)
signature.to_pickle("./output/results/lasso_signature.p")