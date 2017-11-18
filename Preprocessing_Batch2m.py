
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

gse2 = GEOparse.get_GEO(filepath="./data/GEO_data/GSE38337_family.soft.gz")


# In[19]:

plats_2 = list(gse2.gpls.keys())[0]
#print(plats_3)


# Annotation table

# In[20]:
#tenere tessuto, cbmi
#fare dropna (gene symbol e entrez id). nel groupby sul gene symbol e fare la mean
samples2 = gse2.phenotype_data[["characteristics_ch1.5.treatment","characteristics_ch1.6.treatment duration", "characteristics_ch2.4.tissue"]]
samples2 = pd.DataFrame(samples2);
samples2.rename(columns={"characteristics_ch1.5.treatment":"diet","characteristics_ch1.6.treatment duration":"duration", "characteristics_ch2.4.tissue":"tissue"}, inplace=True)
samples2 = samples2.loc[samples2["duration"]=="12 weeks"] 
samples2 = samples2.drop("duration", axis=1)
samples2["cbmi"] = samples2["diet"].apply(lambda x: "obese" if (x == "high-fat diet") else "lean" )
samples2 = samples2.drop("diet", axis=1)


# In[35]:

samples2.to_pickle('./output/results/preprocessing/batch2_m_pheno.p')
with open('./output/results/preprocessing/batch2_m_pheno.txt', 'w') as handle:
    samples2.to_csv(handle, sep='\t')


# # Preprocessing of Expression Data (Batch 4)

# In[23]:
gpl2 = gse2.gpls.get('GPL10333')
samples2_ = gse2.pivot_samples('VALUE')[list(samples2.index)]
#print(samples3_.head())
#samples2 = gse2.pivot_and_annotate("VALUE", gpl=gpl2, annotation_column=["GENE_SYMBOL","GENE"])

# In[24]:
samples2_ann = samples2_.reset_index().merge(gse2.gpls['GPL10333'].table[["ID", "GENE", "GENE_SYMBOL"]],left_on='ID_REF', right_on="ID").set_index('ID_REF')
samples2_ann = samples2_ann.drop("ID", axis=1)
samples2_ann = samples2_ann.dropna()
#del samples3_ann["ID"]
samples2_ann.rename(columns={"GENE":"ENTREZ_GENE_ID"}, inplace=True)

samples2_ann['ENTREZ_GENE_ID'] = samples2_ann['ENTREZ_GENE_ID'].astype(int).astype(str)
#print(samples3_ann.head(), '\n\n', samples3_ann.shape[0])


# # 5) Probes that cannot be mapped to a unique Entrez ID label are excluded from the analysis, as well as those that cannot be mapped to any Entrez ID label at all.

# In[25]:

#print(samples3_ann.shape[0])
samples2_ann = samples2_ann[~samples2_ann.ENTREZ_GENE_ID.str.contains("///")]
samples2_ann = samples2_ann[~samples2_ann.GENE_SYMBOL.str.contains("///")]


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

exprs_2 = samples2_ann.groupby('GENE_SYMBOL').mean()
#print(exprs_3.head())
#print('\n', exprs_3.shape)


# # Write the eprs4 dataframe in a text file

# In[36]:

exprs_2 = exprs_2.T
exprs_2.to_pickle('./output/results/preprocessing/batch2_m_geno.p')
with open('./output/results/preprocessing/batch2_m_geno.txt', 'w') as handle:
    exprs_2.to_csv(handle, sep='\t')

