
# coding: utf-8
    
    # # Batch 1 Preprocessing
    # Group 4: Damiano Chini, Riccardo Gilmozzi, Gianmarco Piccinno & Alessandro Rizzuto
    # Useful links:
    # https://github.com/ComplexityBiosystems/obesity-score
    # https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE2508
    
    # This code operates the preprocessing steps, namely (from the paper):
    #
    # 1) Probes containing missing values are excluded from the analysis.
    #
# 2) Probes are mapped to Entrez ID labels if they are available in the associated platform.
#
# 3) Values corresponding to raw expression counts or gene expression intensity are log2 transformed (if necessary).
#
# 4) Probes mapping to the same Entrez ID label are averaged out.
#
# 5) Probes that cannot be mapped to a unique Entrez ID label are excluded from the analysis, as well as those that cannot be mapped to any Entrez ID label at all.
#
# 6) We apply a simple L 1 normalization in linear space, imposing that the sum of expression of all genes is constant among samples.
#
# After these steps, each data set or batch is represented by a single expression matrix X. Each entry X i j represents the log 2 of the expression intensity of gene i in sample j.
#

# In[2]:

import GEOparse
import pandas as pd
import numpy as np
from functools import *
import re


# In[3]:

gse1 = GEOparse.get_GEO(filepath="./data/GEO_data/GSE2508_family.soft.gz")


# In[4]:

plats_1 = list(gse1.gpls.keys())
#print(plats_1)


# # Phenotype Table (Clinical Data)

# In[5]:

samples1 = gse1.phenotype_data[["platform_id", "title"]]
sample1 = samples1.groupby(["platform_id"]); sample1.groups
d = {}
for l in plats_1:
    #print("\nPlatform: "+str(l)+"\n", sample1.get_group(l))
    #print("\nPlatform: "+str(l)+"\n", sample1.get_group(l)['title'])
    ls = "".join(list(sample1.get_group(l)['title']))
    lf = re.findall("Lean F", ls)
    of = re.findall("Obese F", ls)
    lm = re.findall("Lean M", ls)
    om = re.findall("Obese M", ls)
    #print("LF: ", len(lf), "\nOF: ", len(of), "\nLM: ", len(lm), "\nOM: ", len(om))
    d[l] = {"LF": len(lf), "OF": len(of), "LM": len(lm), "OM": len(om)}
#print(d)
#df = pd.DataFrame(d); print(df)
#df.sum(axis=1)
#df.sum(axis=0)
#print(samples1.head(), len(samples1))


# In[108]:

x = samples1.copy()
x["samples"] = x.index
x["title"] = x['title'].apply(lambda x: x[:-len(x.split()[-1])].strip()).to_frame('samples')
x['gender'] = x['title'].map(lambda x: x.split(' ')[1])
x['cbmi'] = x['title'].map(lambda x: x.split(' ')[0].lower())

# In[20]:

grouped = x.groupby("title")
l = pd.DataFrame.from_dict(grouped.groups)

# In[109]:
print(x)
y = x[["title", "gender", "cbmi"]]

#y = y[["gender", "cbmi"]]
y = y.drop_duplicates("title")
#print(y.shape)
y.index = y.title; y = y[["gender", "cbmi"]]
#print(y); print(y.shape)

y.to_pickle('./output/replication/preprocessing/batch1_pheno.p')
with open('./output/replication/preprocessing/batch1_pheno.txt', 'w') as handle:
    y.to_csv(handle, sep='\t')


# # Preprocessing of Expression Data (Batch 1)

# Batch 1 is composed of five different datasets that use 5 different Affymetrix platforms, each one represents a technical replicate.

# In[9]:

plats_1 = list(gse1.gpls.keys())
#plats_1


# In[10]:

d = {}
samples1 = gse1.phenotype_data[["platform_id", "title"]]
sample1 = samples1.groupby(["platform_id"])

for plat in plats_1:
    d[plat] = gse1.pivot_samples('VALUE')[list(sample1.get_group(plat)[["title"]].index)]

# In[11]:

d_ann = {}
for key in d.keys():
    d_ann[key] = d[key].reset_index().merge(gse1.gpls[key].table[["ID", "ENTREZ_GENE_ID"]],
                                left_on='ID_REF', right_on="ID").set_index('ID_REF')
    d_ann[key].drop(['ID'], axis=1, inplace=True)
    #print(d_ann[key].shape[0])


# # 1) Probes containing missing values are excluded from the analysis.

# In[12]:

for key in d_ann.keys():
    d_ann[key].dropna(inplace=True)
    #print(str(key)+': ', d_ann[key].shape[0])

#print(d_ann['GPL92'].head())


# # 5) Probes that cannot be mapped to a unique Entrez ID label are excluded from the analysis, as well as those that cannot be mapped to any Entrez ID label at all.

# In[13]:

for key in d_ann.keys():
    idx = d_ann[key].ENTREZ_GENE_ID.str.contains("///")#.index
    idx = idx[idx==True].index
    d_ann[key] = d_ann[key].drop(list(idx), axis=0)
    d_ann[key] = d_ann[key].groupby("ENTREZ_GENE_ID").mean()
    #print(d_ann[key].head())
    #print(str(key)+': ', d_ann[key].shape[0])


# # 2) Probes are mapped to Entrez ID labels if they are available in the associated platform.

#
# # 3) Values corresponding to raw expression counts or gene expression intensity are log2 transformed (if necessary).

# In[14]:

for key in d_ann.keys():
    d_ann[key] = np.log2(d_ann[key])
    #print(str(key)+':\n', d_ann[key].head())



# In[61]:

s = [d_ann[key] for key in d_ann.keys()]


# In[38]:

df_final = reduce(lambda left,right: pd.merge(left,right, left_index=True, right_index=True, how='outer'), s)
#print(df_final.head())
#print(df_final.shape)


# # Aggregate data about the 39 samples (tot=195)

# In[22]:

#print(l.head())


# In[66]:

df_c = {}

for column in l.columns:
    df_c[column] = df_final[list(l[column])]

for key in df_c.keys():
    df_c[key] = df_c[key].mean(axis=1)
    df_c[key] = df_c[key].to_frame()
    df_c[key] = df_c[key].rename(columns={0: key})
#print(df_c)

p = [df_c[key] for key in df_c.keys()]

#print(p)
df_f = reduce(lambda left,right: pd.merge(left,right, left_index=True, right_index=True, how='outer'), p)
#print(df_f.shape)
#print(df_f.head())


# # 4) Probes mapping to the same Entrez ID label are averaged out.

# # Write the eprs4 dataframe in a text file

# In[111]:

exprs_1 = df_f.T
exprs_1.to_pickle('./output/replication/preprocessing/batch1_geno.p')
with open('./output/replication/preprocessing/batch1_geno.txt', 'w') as handle:
    exprs_1.to_csv(handle, sep='\t')


# # Comparison with the data provided by the authors

# In[19]:

or_data1 = pd.read_pickle('./data/Paper_data/GSE2508_geno.p')
with open('./output/replication/Comparison.txt', 'w') as handle:
    handle.write("Batch"+"\t"+"#Cols_Pap"+"\t"+"#Rows_Pap"+"\t"+"#Cols_Res"+"\t"+"#Rows_Res"+"\t"+"Rw_Intersection"+"\t"+"Difference"+"\t"+"Jaccard-Diff"\
    "\n"+"Batch1"+"\t"+str(or_data1.shape[0])+"\t"+str(or_data1.shape[1])+"\t"+str(exprs_1.shape[0])\
    +"\t"+str(exprs_1.shape[1])+"\t"+str(len(set.intersection(set(or_data1.columns), set(exprs_1.columns))))\
    +"\t"+str(max(len(set(exprs_1.columns)), len(set(or_data1.columns)))-len(set.intersection(set(exprs_1.columns), set(or_data1.columns))))\
    +"\t"+str(len(set.intersection(set(exprs_1.columns), set(or_data1.columns)))/len(set.union(set(exprs_1.columns), set(or_data1.columns)))))
