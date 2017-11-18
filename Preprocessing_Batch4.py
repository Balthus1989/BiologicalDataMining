
# coding: utf-8

# # Batch 4 Preprocessing
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

gse4 = GEOparse.get_GEO(filepath="./data/GEO_data/GSE48964_family.soft.gz")


# In[3]:

plats_4 = list(gse4.gpls.keys())[0]
#print(plats_4)


# Annotation table

# In[20]:

samples4 = gse4.phenotype_data["source_name_ch1"]
samples4 = pd.DataFrame(samples4); samples4.head()
samples4.rename(columns={"source_name_ch1":"cbmi"}, inplace=True); samples4.head()
samples4["cbmi"] = samples4["cbmi"].apply(lambda x: 'obese' if x.split(' ')[1].lower()=='obese' else 'lean'); samples4.head()
#print(samples4.head()); print(len(samples4))
#print(samples4.shape)


# In[25]:

samples4.to_pickle('./output/replication/preprocessing/batch4_pheno.p')
with open('./output/replication/preprocessing/batch4_pheno.txt', 'w') as handle:
    samples4.to_csv(handle, sep='\t')


# # Preprocessing of Expression Data (Batch 4)

# In[6]:

samples4_exprs = gse4.pivot_samples('VALUE')[list(samples4.index)]
samples4_exprs.index.name = 'ID'
#print('Expression Table', '\n', samples4_exprs.head())


# We are considering the annotation table the one contained in the GPL6244.annot file (https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GPL6244).
#
# ID = ID from Platform data table
#
# Gene ID = Entrez Gene identifier
#
# The column containing the Entrez_ID is Gene ID!!!

# In[7]:

ann_table4 = pd.read_csv('./data/GEO_data/GPL6244.annot', sep='\t', skiprows=27, dtype=str, na_values='NaN',usecols=[0,3], index_col=0)#.dropna()
#print('Annotation Table: ', '\n', ann_table4.head())


# Remove probes without ENTREZ ID

# # 1) Probes containing missing values are excluded from the analysis.

# In[8]:

#print(samples4_exprs.shape)
samples4_exprs = samples4_exprs.dropna()
#print(samples4_exprs.shape)


# # 5) Probes that cannot be mapped to a unique Entrez ID label are excluded from the analysis, as well as those that cannot be mapped to any Entrez ID label at all.

# In[9]:

ann_table = ann_table4.dropna()
#print(ann_table.head())


# Remove ids that refer to different Entrez ids:
# for example
#          ID                              Gene ID
# 2   7896740                81099///79501///26682
#

# In[10]:

ann_table = ann_table[~ann_table['Gene ID'].str.contains("///")]
#print(ann_table.head())
#print(ann_table.shape)


# # 2) Probes are mapped to Entrez ID labels if they are available in the associated platform.

# In[11]:

#print(type(ann_table.index[0]))
#print(type(samples4_exprs.index[0]))
samples4_exprs.index = samples4_exprs.index.astype(str)
#print(type(ann_table.index[0]))
#print(type(samples4_exprs.index[0]))


# In[12]:

exprs_4 = ann_table.merge(samples4_exprs, left_index=True, right_index=True, how='inner')
#exprs_4.index = exprs_4['Gene ID']; del exprs_4['Gene ID']
#print(exprs_4.head())
#print('\nShape of the complete DataFrame: ', exprs_4.shape)


#
# # 3) Values corresponding to raw expression counts or gene expression intensity are log2 transformed (if necessary).

# # 4) Probes mapping to the same Entrez ID label are averaged out.

# In[13]:

exprs_4 = exprs_4.groupby('Gene ID').mean()
#print(exprs_4.head())
#print('\n', exprs_4.shape)


# # Write the eprs4 dataframe in a text file

# In[24]:

exprs_4 = exprs_4.T
exprs_4.to_pickle('./output/replication/preprocessing/batch4_geno.p')
with open('./output/replication/preprocessing/batch4_geno.txt', 'w') as handle:
    exprs_4.to_csv(handle, sep='\t')


# # Comparison with the data provided by the authors

# In[18]:

or_data4 = pd.read_pickle('./data/Paper_data/GSE48964_geno.p')

with open('./output/replication/Comparison.txt', 'a') as handle:
    handle.writelines("\nBatch4"+"\t"+str(or_data4.shape[0])+"\t"+str(or_data4.shape[1])+"\t"+str(exprs_4.shape[0])\
    +"\t"+str(exprs_4.shape[1])+"\t"+str(len(set.intersection(set(exprs_4.columns), set(or_data4.columns))))+"\t"\
    +"\t"+str(max(len(set(exprs_4.columns)), len(set(or_data4.columns)))-len(set.intersection(set(exprs_4.columns), set(or_data4.columns))))\
    +"\t"+str(len(set.intersection(set(exprs_4.columns), set(or_data4.columns)))/len(set.union(set(exprs_4.columns), set(or_data4.columns)))))
