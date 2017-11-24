import pandas as pd

list_batch = ['batch1', 'batch2','batch3','batch4']
dictionary = {}
for batch in list_batch:
    dictionary[batch] = pd.read_pickle("./output/results/preprocessing/" + batch + "_geno.p").T

dictionary["batch4"]


hgcn = pd.read_csv("hgnc.txt", sep="\t")
#print(hgcn.head())

hg = hgcn[["Approved Symbol", "Entrez Gene ID(supplied by NCBI)"]]
hg.rename(columns = {"Approved Symbol":"Gene_Symbol", "Entrez Gene ID(supplied by NCBI)":"ENTREZ_ID"}, inplace=True)
hg = hg.dropna()
hg["ENTREZ_ID"] = hg["ENTREZ_ID"].astype(int).astype(str)

for h in dictionary.keys():    
    merge = hg.merge(dictionary[h], left_on="ENTREZ_ID", right_index=True)
    print(merge.head())
    #print("Original shape: "+str(h.shape[0]))
    #print("Definit_shape: "+str(merge.shape[0]))
    
    cols = list(merge.columns)
    cols.remove("ENTREZ_ID")
    merge = merge[cols]
    print(merge.head)
    
    merge.to_csv("./output/results/preprocessing/"+str(h)+"GS_geno.txt", sep="\t", index=False)
    merge.to_pickle("./output/results/preprocessing/"+str(h)+"GS_geno.p")