
# coding: utf-8

# # Batch 6 Preprocessing
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


import GEOparse
import pandas as pd

gse6 = GEOparse.get_GEO(filepath="./data/GEO_data/GSE64567_family.soft.gz")

plats_6 = list(gse6.gpls.keys())[0]
#print(plats_6)


# Annotation table

samples6 = gse6.phenotype_data["characteristics_ch1.4.bmi (kg/m2)"]
samples6 = pd.DataFrame(samples6); samples6.head()
samples6.rename(columns={"characteristics_ch1.4.bmi (kg/m2)":"bmi"}, inplace=True)

samples6["cbmi"] = samples6["bmi"].apply(lambda x: "obese" if (float(x) > 30) else ("lean" if (float(x) <= 25) else ("overweight" if (float(x) > 25) & (float(x) <= 30) else "STRANGE")) )

samples6 = samples6[['cbmi', 'bmi']]
#print(samples6)

samples6.to_pickle('./output/results/preprocessing/batch6_pheno.p')
with open('./output/results/preprocessing/batch6_pheno.txt', 'w') as handle:
    samples6.to_csv(handle, sep='\t')


# # Preprocessing of Expression Data (Batch 6)

samples6_ = gse6.pivot_samples('VALUE')[list(samples6.index)]
#print(samples6_.head())

samples6_ann = samples6_.reset_index().merge(gse6.gpls['GPL10558'].table[["ID", "Symbol", "Entrez_Gene_ID"]],
                                left_on='ID_REF', right_on="ID").set_index('ID_REF')
#print(samples6_ann.head())

del samples6_ann["ID"]
samples6_ann['Entrez_Gene_ID'] = samples6_ann['Entrez_Gene_ID'].astype(str)
#print(samples6_ann.head(), '\n\n', samples6_ann.shape[0])


# # 5) Probes that cannot be mapped to a unique Entrez ID label are excluded from the analysis, as well as those that cannot be mapped to any Entrez ID label at all.

samples6_ann = samples6_ann[~samples6_ann.Entrez_Gene_ID.str.contains("///")]
#print(samples6_ann.shape[0])

samples6_ann['Entrez_Gene_ID'] = samples6_ann['Entrez_Gene_ID'].astype(float)
samples6_ann = samples6_ann.dropna()
#print(samples6_ann.shape[0])

samples6_ann['Entrez_Gene_ID'] = samples6_ann['Entrez_Gene_ID'].astype(int)
samples6_ann['Entrez_Gene_ID'] = samples6_ann['Entrez_Gene_ID'].astype(str)


# # 1) Probes containing missing values are excluded from the analysis.
# # 2) Probes are mapped to Entrez ID labels if they are available in the associated platform.

#
# # 3) Values corresponding to raw expression counts or gene expression intensity are log2 transformed (if necessary).

# # 4) Probes mapping to the same Entrez ID label are averaged out.

exprs_6 = samples6_ann.groupby('Symbol').mean()
#print(exprs_6.head())
#print('\n', exprs_6.shape)


# # Write the eprs6 dataframe in a text file

exprs_6 = exprs_6.T
exprs_6.to_pickle('./output/results/preprocessing/batch6_geno.p')
with open('./output/results/preprocessing/batch6_geno.txt', 'w') as handle:
    exprs_6.to_csv(handle, sep='\t')
