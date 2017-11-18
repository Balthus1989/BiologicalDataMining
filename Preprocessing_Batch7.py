
# coding: utf-8

# # Batch 3 Preprocessing
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

# In[17]:

import GEOparse
import pandas as pd


# In[18]:

gse7 = GEOparse.get_GEO(filepath="./data/GEO_data/GSE65540_family.soft.gz")


# In[19]:

plats_7 = list(gse7.gpls.keys())[0]
#print(plats_3)


# Annotation table

# In[20]:

samples3 = gse3.phenotype_data["characteristics_ch1.1.bmi"]
samples3 = pd.DataFrame(samples3); samples3.head()
samples3.rename(columns={"characteristics_ch1.1.bmi":"bmi"}, inplace=True)
samples3["cbmi"] = samples3["bmi"].apply(lambda x: "obese" if (float(x) > 30) else ("lean" if (float(x) <= 25) else ("overweight" if (float(x) > 25) & (float(x) <= 30) else "STRANGE")) )
samples3 = samples3[['cbmi', 'bmi']]

#print(samples3)


# In[21]:

or_samples3 = pd.read_pickle('./data/Paper_data/GSE27949_pheno.p')
#print(type(or_samples3['bmi'][0]), type(or_samples3['cbmi'][0]), type(or_samples3.index[0]))
#print(type(samples3['bmi'][0]), type(samples3['cbmi'][0]), type(samples3.index[0]))
or_samples3['bmi'] = or_samples3['bmi'].astype(str)
#print(type(or_samples3['bmi'][0]), type(or_samples3['cbmi'][0]), type(or_samples3.index[0]))
#print(type(samples3['bmi'][0]), type(samples3['cbmi'][0]), type(samples3.index[0]))

or_samples3_d = or_samples3.to_dict()
samples3_d = samples3.to_dict()
or_samples3_d == samples3_d


# In[35]:

samples3.to_pickle('./output/replication/preprocessing/batch3_pheno.p')
with open('./output/replication/preprocessing/batch3_pheno.txt', 'w') as handle:
    samples3.to_csv(handle, sep='\t')


# # Preprocessing of Expression Data (Batch 4)

# In[23]:

samples3_ = gse3.pivot_samples('VALUE')[list(samples3.index)]
#print(samples3_.head())


# In[24]:

samples3_ann = samples3_.reset_index().merge(gse3.gpls['GPL570'].table[["ID", "ENTREZ_GENE_ID"]],
                                left_on='ID_REF', right_on="ID").set_index('ID_REF')
del samples3_ann["ID"]
samples3_ann['ENTREZ_GENE_ID'] = samples3_ann['ENTREZ_GENE_ID'].astype(str)
#print(samples3_ann.head(), '\n\n', samples3_ann.shape[0])


# # 5) Probes that cannot be mapped to a unique Entrez ID label are excluded from the analysis, as well as those that cannot be mapped to any Entrez ID label at all.

# In[25]:

#print(samples3_ann.shape[0])
samples3_ann = samples3_ann[~samples3_ann.ENTREZ_GENE_ID.str.contains("///")]


# In[26]:

samples3_ann['ENTREZ_GENE_ID'] = samples3_ann['ENTREZ_GENE_ID'].astype(float)
#print(samples3_ann.shape[0])
samples3_ann = samples3_ann.dropna()
#print(samples3_ann.shape[0])
samples3_ann['ENTREZ_GENE_ID'] = samples3_ann['ENTREZ_GENE_ID'].astype(int)
samples3_ann['ENTREZ_GENE_ID'] = samples3_ann['ENTREZ_GENE_ID'].astype(str)


# # 1) Probes containing missing values are excluded from the analysis.

# In[27]:

#print(samples3_ann.shape)
samples3_ann = samples3_ann.dropna()
#print(samples3_ann.shape)


# # 2) Probes are mapped to Entrez ID labels if they are available in the associated platform.

#
# # 3) Values corresponding to raw expression counts or gene expression intensity are log2 transformed (if necessary).

# # 4) Probes mapping to the same Entrez ID label are averaged out.

# In[28]:

exprs_3 = samples3_ann.groupby('ENTREZ_GENE_ID').median()
#print(exprs_3.head())
#print('\n', exprs_3.shape)


# # Write the eprs4 dataframe in a text file

# In[36]:

exprs_3 = exprs_3.T
exprs_3.to_pickle('./output/replication/preprocessing/batch3_geno.p')
with open('./output/replication/preprocessing/batch3_geno.txt', 'w') as handle:
    exprs_3.to_csv(handle, sep='\t')


# # Comparison with the data provided by the authors

# In[32]:

or_data3 = pd.read_pickle('./data/Paper_data/GSE27949_geno.p')

with open('./output/replication/Comparison.txt', 'a') as handle:
    handle.writelines("\nBatch3"+"\t"+str(or_data3.shape[0])+"\t"+str(or_data3.shape[1])+"\t"+str(exprs_3.shape[0])\
    +"\t"+str(exprs_3.shape[1])+"\t"+str(len(set.intersection(set(or_data3.columns), set(exprs_3.columns))))\
    +"\t"+str(max(len(set(exprs_3.columns)), len(set(or_data3.columns)))-len(set.intersection(set(exprs_3.columns), set(or_data3.columns))))\
    +"\t"+str(len(set.intersection(set(exprs_3.columns), set(or_data3.columns)))/len(set.union(set(exprs_3.columns), set(or_data3.columns)))))
    
