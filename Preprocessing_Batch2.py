
# coding: utf-8

# # Batch 2 Preprocessing
# Group 4: Damiano Chini, Riccardo Gilmozzi, Gianmarco Piccinno & Alessandro Rizzuto
# Useful links:
# https://github.com/ComplexityBiosystems/obesity-score
# https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE48964

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

# In[1]:

import GEOparse
import pandas as pd


# In[2]:

gse2 = GEOparse.get_GEO(filepath="./data/GEO_data/GSE26637_family.soft.gz")


# In[3]:

plats_2 = list(gse2.gpls.keys())[0]
#print(plats_2)


# Annotation table

# In[4]:

samples2 = gse2.phenotype_data[["characteristics_ch1.0.gender", "characteristics_ch1.2.stimulation", "characteristics_ch1.3.resistance status"]]
samples2 = samples2.rename(columns={'characteristics_ch1.0.gender':'gender', 'characteristics_ch1.2.stimulation':'fasting_status',
                         'characteristics_ch1.3.resistance status':'insulin_status'})
samples2['cbmi'] = samples2['insulin_status'].map(lambda x: 'lean' if x == 'sensitive' else 'obese')
#print(samples2)


# In[15]:

samples2.to_pickle('./output/replication/preprocessing/batch2_pheno.p')
with open('./output/replication/preprocessing/batch2_pheno.txt', 'w') as handle:
    samples2.to_csv(handle, sep='\t')


# # Preprocessing of Expression Data (Batch 2)

# In[6]:

samples2_exprs = gse2.pivot_samples('VALUE')[list(samples2.index)]
#print('Expression Table', '\n', samples2_exprs.head())


# In[7]:

samples2_ann = samples2_exprs.reset_index().merge(gse2.gpls['GPL570'].table[["ID", "ENTREZ_GENE_ID"]],
                                left_on='ID_REF', right_on="ID").set_index('ID_REF')
samples2_ann.drop('ID', inplace=True, axis=1)
samples2_ann['ENTREZ_GENE_ID'] = samples2_ann['ENTREZ_GENE_ID'].astype(str)
#print(samples2_ann.head())


# # 5) Probes that cannot be mapped to a unique Entrez ID label are excluded from the analysis, as well as those that cannot be mapped to any Entrez ID label at all.

# In[8]:

#print(samples2_ann.shape[0])
samples2_ann = samples2_ann[~samples2_ann.ENTREZ_GENE_ID.str.contains("///")].dropna()
samples2_ann['ENTREZ_GENE_ID'].astype(float, inplace=True)
#print(samples2_ann.shape[0])
samples2_ann = samples2_ann.dropna()
#print(samples2_ann.shape[0])


# # 1) Probes containing missing values are excluded from the analysis.

# # 2) Probes are mapped to Entrez ID labels if they are available in the associated platform.

#
# # 3) Values corresponding to raw expression counts or gene expression intensity are log2 transformed (if necessary).

# # 4) Probes mapping to the same Entrez ID label are averaged out.

# In[9]:

exs_2 = samples2_ann.groupby('ENTREZ_GENE_ID').median()
#print(exprs_2.head())
#print('\n', exprs_2.shape)


# # Write the eprs4 dataframe in a text file

# In[16]:

exprs_2 = exs_2.T
exprs_2.to_pickle('./output/replication/preprocessing/batch2_geno.p')
with open('./output/replication/preprocessing/batch2_geno.txt', 'w') as handle:
    exprs_2.to_csv(handle, sep='\t')


# # Comparison with the data provided by the authors

# In[12]:

or_data2 = pd.read_pickle('./data/Paper_data/GSE26637_geno.p')

with open('./output/replication/Comparison.txt', 'a') as handle:
    handle.writelines("\nBatch2"+"\t"+str(or_data2.shape[0])+"\t"+str(or_data2.shape[1])+"\t"+str(exprs_2.shape[0])\
    +"\t"+str(exprs_2.shape[1])+"\t"+str(len(set.intersection(set(or_data2.columns), set(exprs_2.columns))))\
    +"\t"+str(max(len(set(exprs_2.columns)), len(set(or_data2.columns)))-len(set.intersection(set(exprs_2.columns), set(or_data2.columns))))
    +"\t"+str(len(set.intersection(set(exprs_2.columns), set(or_data2.columns)))/len(set.union(set(exprs_2.columns), set(or_data2.columns)))))
