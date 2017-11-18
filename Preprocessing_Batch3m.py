import GEOparse
import pandas as pd
import os

os.chdir("/home/gianmarco/Scrivania/First_Semester/Data_Mining_Laboratory/Project")

gse3m = GEOparse.get_GEO(filepath="./new_data/GSE39549_family.soft.gz")

plats_3m = list(gse3m.gpls.keys())[0]
#print(plats_4)


# Annotation table

d = pd.DataFrame(gse3m.phenotype_data[["characteristics_ch1.1.treatment protocol", "characteristics_ch1.2.time"]]).dropna()
d.rename(columns={"characteristics_ch1.1.treatment protocol": "diet", "characteristics_ch1.2.time": "time"}, inplace=True)

d["time"] = d["time"].map(lambda x: x[:-5]).astype(int)
d = d[d["time"]>=8]

anns_3m = {}

for el in d.index:
    if d["diet"][el] == "Normal diet":
        anns_3m[el] = "lean"
    else:
        anns_3m[el] = "obese"

anns_3m = pd.DataFrame([anns_3m])


anns_3m.T.to_pickle('./output/results/preprocessing/batch3m_pheno.p')
with open('./output/results/preprocessing/batch3m_pheno.txt', 'w') as handle:
    anns_3m.T.to_csv(handle, sep='\t', header=False)

#Expression data analysis

samples3m_ = gse3m.pivot_samples('VALUE')[list(anns_3m.columns)]


samples3m_ann = samples3m_.reset_index().merge(gse3m.gpls[plats_3m].table[["ID", "Entrez_Gene_ID", "Symbol"]],
                                left_on = "ID_REF", right_on="ID").set_index('ID').dropna()

exprs_3m = samples3m_ann.groupby('Symbol').mean()
exprs_3m.drop("Entrez_Gene_ID", axis=1, inplace=True)

del exprs_3m.index.name

#exprs_5.boxplot()

exprs_3m.to_pickle('./output/results/preprocessing/batch3m_geno.p')
with open('./output/results/preprocessing/batch3m_geno.txt', 'w') as handle:
    exprs_3m.to_csv(handle, sep='\t')   