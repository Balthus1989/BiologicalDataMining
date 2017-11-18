import GEOparse
import pandas as pd
import os

os.chdir("/home/gianmarco/Scrivania/First_Semester/Data_Mining_Laboratory/Project")

gse5 = GEOparse.get_GEO(filepath="./new_data/GSE62117_family.soft.gz")

plats_5 = list(gse5.gpls.keys())[0]
#print(plats_4)


# Annotation table

d = pd.DataFrame(gse5.phenotype_data["title"])
d["keep"] = gse5.phenotype_data["title"].map(lambda x: x.split("-")[1] =="CNTRL")
ix_s = d[d["keep"] == True]

samples5 = gse5.phenotype_data.loc[list(ix_s.index)]
samples5 = samples5[["characteristics_ch1.3.gender", "characteristics_ch1.5.bmi", ]]
samples5.rename(columns={"characteristics_ch1.3.gender":"gender", "characteristics_ch1.5.bmi":"bmi"}, inplace=True)
samples5["cbmi"] = samples5["bmi"].apply(lambda x: "obese" if (float(x) > 30) else ("lean" if (float(x) <= 25) else ("overweight" if (float(x) > 25) & (float(x) <= 30) else "STRANGE")) )

list(samples5["cbmi"]).count("overweight")
list(samples5["cbmi"]).count("lean")
list(samples5["cbmi"]).count("obese")

samples5 = samples5[~(samples5["cbmi"] == "overweight")]

samples5.to_pickle('./output/results/preprocessing/batch5_pheno.p')
with open('./output/results/preprocessing/batch5_pheno.txt', 'w') as handle:
    samples5.to_csv(handle, sep='\t')


#Expression data analysis

samples5_ = gse5.pivot_samples('VALUE')[list(samples5.index)]


samples5_ann = samples5_.reset_index().merge(gse5.gpls[plats_5].table[["GENE", "GENE_SYMBOL"]],
                                left_on='ID_REF', right_index=True).set_index('ID_REF').dropna()

exprs_5 = samples5_ann.groupby('GENE_SYMBOL').mean()
exprs_5.drop("GENE", axis=1, inplace=True)
del exprs_5.index.name

#exprs_5.boxplot()

exprs_5.to_pickle('./output/results/preprocessing/batch5_geno.p')
with open('./output/results/preprocessing/batch5_geno.txt', 'w') as handle:
    exprs_5.to_csv(handle, sep='\t')