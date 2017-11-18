
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

gse4 = GEOparse.get_GEO(filepath="./data/GEO_data/GSE32095_family.soft.gz")


# In[19]:

plats_4 = list(gse4.gpls.keys())[0]
#print(plats_3)


# Annotation table

# In[20]:
#tenere tessuto, cbmi
#fare dropna (gene symbol e entrez id). nel groupby sul gene symbol e fare la mean
samples4 = gse4.phenotype_data[["characteristics_ch1.0.tissue","characteristics_ch1.2.treatment", "characteristics_ch1.1.genotype"]]
samples4 = pd.DataFrame(samples4);
samples4.rename(columns={"characteristics_ch1.0.tissue":"tissue","characteristics_ch1.2.treatment":"diet", "characteristics_ch1.1.genotype":"mouse_type"}, inplace=True)
samples4 = samples4.loc[samples4["mouse_type"] == "wild-type"]
samples4 = samples4.loc[samples4["tissue"] == "Epididymal adipose"]

samples4 = samples4.drop("mouse_type", axis=1)
samples4["cbmi"] = samples4["diet"].apply(lambda x: "obese" if (x == "high fat diet HFD") else "lean" )
samples4 = samples4.drop("diet", axis=1)


# In[35]:

samples4.to_pickle('./output/results/preprocessing/batch4_m_pheno.p')
with open('./output/results/preprocessing/batch4_m_pheno.txt', 'w') as handle:
    samples4.to_csv(handle, sep='\t')


# # Preprocessing of Expression Data (Batch 4)

# In[23]:
samples4_ = gse4.pivot_samples('VALUE')[list(samples4.index)]
#print(samples3_.head())
#samples4 = gse4.pivot_and_annotate("VALUE", gpl=gpl4, annotation_column=["GENE_SYMBOL","GENE"])

# In[24]:

samples4_ann = gse4.gpls['GPL1261'].table
samples4_ann = samples4_.reset_index().merge(gse4.gpls['GPL1261'].table[["ID", "ENTREZ_GENE_ID", "Gene Symbol"]],left_on='ID_REF', right_on="ID").set_index('ID_REF')
samples4_ann = samples4_ann.drop("ID", axis=1)
samples4_ann = samples4_ann.dropna()
#del samples3_ann["ID"]
#samples4_ann.rename(columns={"GENE":"ENTREZ_GENE_ID"}, inplace=True)

samples4_ann['ENTREZ_GENE_ID'] = samples4_ann['ENTREZ_GENE_ID'].astype(str)
#print(samples3_ann.head(), '\n\n', samples3_ann.shape[0])


# # 5) Probes that cannot be mapped to a unique Entrez ID label are excluded from the analysis, as well as those that cannot be mapped to any Entrez ID label at all.

# In[25]:

#print(samples3_ann.shape[0])
samples4_ann = samples4_ann[~samples4_ann.ENTREZ_GENE_ID.str.contains("///")]
samples4_ann = samples4_ann[~samples4_ann["Gene Symbol"].str.contains("///")]


# In[26]:
'''
samples3_ann['ENTREZ_GENE_ID'] = samples3_ann['ENTREZ_GENE_ID'].astype(float)
#print(samples3_ann.shape[0])
samples3_ann = samples3_ann.dropna()
#print(samples3_ann.shape[0])
samples3_ann['ENTREZ_GENE_ID'] = samples3_ann['ENTREZ_GENE_ID'].astype(int)
samples3_ann['ENTREZ_GENE_ID'] = samples3_ann['ENTREZ_GENE_ID'].astype(str)
'''

# # 1) Probes containing missing values are excluded from the analysis.

# In[27]:
'''
#print(samples3_ann.shape)
samples3_ann = samples3_ann.dropna()
#print(samples3_ann.shape)
'''

# # 2) Probes are mapped to Entrez ID labels if they are available in the associated platform.

#
# # 3) Values corresponding to raw expression counts or gene expression intensity are log2 transformed (if necessary).

# # 4) Probes mapping to the same Entrez ID label are averaged out.

# In[28]:

exprs_4 = samples4_ann.groupby('Gene Symbol').mean()
#print(exprs_3.head())
#print('\n', exprs_3.shape)


# # Write the eprs4 dataframe in a text file

# In[36]:

exprs_4 = exprs_4.T
exprs_4.to_pickle('./output/results/preprocessing/batch4_m_geno.p')
with open('./output/results/preprocessing/batch4_m_geno.txt', 'w') as handle:
    exprs_4.to_csv(handle, sep='\t')

