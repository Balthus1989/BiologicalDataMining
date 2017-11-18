import numpy as np
import pandas as pd
from sklearn.decomposition import PCA

import matplotlib.pyplot as plt
#get_ipython().magic('matplotlib inline')
import seaborn as sns

from scipy.stats import ks_2samp

import copy

import os
os.chdir("/home/gianmarco/Scrivania/First_Semester/Data_Mining_Laboratory/Project")

human_batches = ["batch1", "batch2", "batch3", "batch4", "batch5", "batch6"]




def normalize_df(df):
    """
    input: pandas dataframe, values are log2 gene expression, genes in columns, samples in rows.
    output: same, with sum of genes equal 1 for all samples (in linear space).
    """

    # exponentiate
    x = pow(2.,df)

    # rows sum to 1
    x = x.div(x.sum(axis=1), axis=0)

    return np.log2(x)



dfs = {}
cols = {}
for batch in human_batches:
    if batch != "batch6":
        dfs[batch] = pd.read_csv("./output/results/preprocessing/"+batch+"_geno.txt", sep="\t", index_col=0)
        #cols[batch] = {}
        #for col in dfs[batch].columns:     
        #    cols[batch][col] = float
    else:
        dfs[batch] = pd.read_csv("./output/results/preprocessing/"+batch+"_geno.txt", sep="\t", index_col=0, header="infer").T
        #cols[batch] = {}
        #for col in dfs[batch].columns:     
        #    cols[batch][col] = float    
           
for batch in dfs:           
    dfs[batch]=(normalize_df(dfs[batch]))
    
geno = pd.concat(dfs, axis=1).dropna(axis=0).T
geno.index = geno.index.droplevel(0)

#

dfs = {}
for gse in human_batches:
    tmp = pd.read_pickle("./output/results/preprocessing/"+gse+"_pheno.p")
    tmp["batch"] = gse
    dfs[gse]=(tmp)
pheno = pd.concat(dfs)
pheno.index = pheno.index.droplevel(0)

geno.drop(geno.columns[geno.isnull().sum()!=0],axis=1,inplace=True)

print("batch \t lean \t overw \t obese \t total")
print("---")
for gse in human_batches:
    l = ((pheno["batch"]==gse) & (pheno["cbmi"]=="lean")).sum()
    over = ((pheno["batch"]==gse) & (pheno["cbmi"]=="overweight")).sum()
    o = ((pheno["batch"]==gse) & (pheno["cbmi"]=="obese")).sum()
    print("%s \t %d \t %d \t %d \t %d" % (gse,l,over,o,l+over+o))
l =  (pheno["cbmi"]=="lean").sum()
over = ((pheno["cbmi"]=="overweight")).sum()
o = ((pheno["cbmi"]=="obese")).sum()
print("---")
print("%s \t %d \t %d \t %d \t %d" % ("ALL\t",l,over,o,l+over+o))

idx = pheno[(pheno["cbmi"]=="overweight")].index

geno.drop(idx,axis=0,inplace=True)
pheno.drop(idx,axis=0,inplace=True)

print("genes:",geno.shape[1])
print("samples (lean):",geno.loc[pheno.cbmi=="lean"].shape[0])
print("samples (obese):",geno.loc[pheno.cbmi=="obese"].shape[0])

ogeno = copy.deepcopy(geno)

def find_level(df,pheno,n=10):
    """
    returns 1./pvalue of distro of n first components being from different distros when split by cbmi
    args:
    df = dataframe with expression values
    pheno = dataframe with labels and metadata
    n = max numb of components
    """
    # make sure we have enough dimensions
    n=min(n,df.shape[0])

    res = []
    pca = PCA()
    trans = pca.fit_transform(df)
    df_trans = pd.DataFrame(data = trans[:,:n],index = df.index)

    for i in range(n):
        res.append([i,1/ks_2samp(
            df_trans.loc[(pheno.cbmi=="lean") | (pheno.cbmi=="overweight")].iloc[:,i],
            df_trans.loc[pheno.cbmi=="obese"].iloc[:,i]
        ).pvalue])
    return np.array(res).T

plt.figure(figsize=(17,3))
for i,batch in enumerate(human_batches):
    plt.subplot(1,6,i+1)
    plt.title(batch,size=14)
    plt.xlabel("eigengene index",size=14)
    if i==0: plt.ylabel("$p$-value, KS lean/obese",size=14)
    y = np.log10(find_level(geno.loc[pheno.batch==batch],pheno)[1])
    plt.plot(range(1,y.shape[0]+1),y)
    yticks = plt.yticks()[0][::2]
    plt.yticks(yticks,["$10^{-%1.1f}$"%i for i in yticks],fontsize=14)
    plt.subplots_adjust(wspace=0.3)
plt.savefig("./output/replication/our/figures/FigureP1.pdf")

pca = PCA()

for batch in human_batches:
    effect_strength = find_level(geno.loc[pheno.batch==batch],pheno)[1]
    while np.argmax(effect_strength)!=0:
        tmp = pca.fit_transform(geno.loc[pheno.batch==batch])
        tmp[:,0]=0
        print("One principal component was set to zero on batch",batch)
        geno.loc[pheno.batch==batch] = pca.inverse_transform(tmp)
        effect_strength = find_level(geno.loc[pheno.batch==batch],pheno,n=min(10,geno.loc[pheno.batch==batch].shape[0]))[1]

plt.figure(figsize=(17,3))
for i,batch in enumerate(human_batches):
    plt.subplot(1,6,i+1)
    plt.title(batch,size=14)
    plt.xlabel("eigengene index",size=14)
    if i==0: plt.ylabel("$p$-value, KS lean/obese",size=14)
    y = np.log10(find_level(geno.loc[pheno.batch==batch],pheno)[1])
    plt.plot(range(1,y.shape[0]+1),y)
    yticks = plt.yticks()[0][::2]
    plt.yticks(yticks,["$10^{-%1.1f}$"%i for i in yticks],fontsize=14)
    plt.subplots_adjust(wspace=0.3)
plt.savefig("./output/replication/our/figures/FigureP2.pdf")

plt.title("ALL GSE's",size=14)
plt.plot(range(1,11),np.log10(find_level(geno,pheno)[1]))
plt.xlabel("eigegene index",size=14)
plt.ylabel("$p$-value (KS lean/obese)",size=14)
plt.xticks(fontsize=14)
yticks = plt.yticks()[0][::2]
plt.yticks(yticks,["$10^{-%d}$"%i for i in yticks],fontsize=16)
#plt.show()
eigvect_k = np.argmax(find_level(geno,pheno)[1])
print("Picking eigengene number",eigvect_k+1)

pca = PCA()
trans = pca.fit_transform(geno)
trans[:,:eigvect_k]=0
n_geno = pd.DataFrame(index = geno.index,columns=geno.columns,data=pca.inverse_transform(trans))

K=7

pca_geno = pd.DataFrame(
    data = PCA(n_components=K,whiten=True).fit_transform(geno),
    index = geno.index).T.corr()

pca_ogeno = pd.DataFrame(
    data = PCA(n_components=K,whiten=True).fit_transform(ogeno),
    index = ogeno.index).T.corr()

pca_ngeno = pd.DataFrame(
    data = PCA(n_components=K,whiten=True).fit_transform(n_geno),
    index = ogeno.index).T.corr()

matrix_colors_dict = dict(zip(np.unique(pheno.batch),
                              np.array([sns.palettes.light_palette(x)[2] for x in sns.color_palette(n_colors=9)])[[3,4,5,6,7, 8]]
                             ))
matrix_colors_dict["lean"] = sns.palettes.light_palette("green")[1]
matrix_colors_dict["obese"] = sns.palettes.light_palette("red")[1]

# sort by batch
order = pheno.sort_values(by=["batch","cbmi"][::1]).index
batchcolors_bybatch = [matrix_colors_dict[pheno.loc[x,"batch"]] for x in order]
cbmicolors_bybatch = [matrix_colors_dict[pheno.loc[x,"cbmi"]] for x in order]

# sort by cbmi
order = pheno.sort_values(by=["batch","cbmi"][::-1]).index
batchcolors_bycbmi = [matrix_colors_dict[pheno.loc[x,"batch"]] for x in order]
cbmicolors_bycbmi = [matrix_colors_dict[pheno.loc[x,"cbmi"]] for x in order]

batch_sizes = np.array([39,20,23,6])
cbmi_size = np.array([39,49])


order = pheno.sort_values(by=["batch","cbmi"][::1]).index
big_ax = sns.clustermap(pca_geno.loc[order,order],
               row_cluster=False,col_cluster=False,
               row_colors=cbmicolors_bybatch,
               col_colors=batchcolors_bybatch,
               xticklabels=False,
               yticklabels=False,
               figsize=(6,6)
           )
ax = big_ax.ax_heatmap
big_ax.cax.set_visible(False)
N = n_geno.shape[(0)]
nn=0
for b,n in enumerate(batch_sizes[:-1]):
    nn=nn+n
    ax.axhline(N-nn,color="black")
    ax.axvline(nn,color="black")
    ax.text(nn-n/2,N+1,"Batch"+str(b+1),rotation=0,fontsize=13.5,color="black",horizontalalignment="center")
ax.text(N-3,N+1,"B4",rotation=0,fontsize=13.5,color="black",horizontalalignment="center")
#plt.savefig("./output/replication/our/figures/Figure1a.pdf")
plt.show()


order = pheno.sort_values(by=["batch","cbmi"][::1]).index
big_ax = sns.clustermap(pca_ngeno.loc[order,order],
               row_cluster=False,col_cluster=False,
               row_colors=cbmicolors_bybatch,
               col_colors=batchcolors_bybatch,
               xticklabels=False,
               yticklabels=False,
               figsize=(6,6)
           )
ax = big_ax.ax_heatmap
big_ax.cax.set_visible(False)
N = n_geno.shape[(0)]
nn=0
for b,n in enumerate(batch_sizes[:-1]):
    nn=nn+n
    ax.axhline(N-nn,color="black")
    ax.axvline(nn,color="black")
    ax.text(nn-n/2,N+1,"Batch"+str(b+1),rotation=0,fontsize=13.5,color="black",horizontalalignment="center")
ax.text(N-3,N+1,"B4",rotation=0,fontsize=13.5,color="black",horizontalalignment="center")
#plt.savefig("./output/replication/our/figures/Figure1b.pdf")
plt.show()



order = pheno.sort_values(by=["batch","cbmi"][::-1]).index
big_ax = sns.clustermap(pca_ngeno.loc[order,order],
               row_cluster=False,col_cluster=False,
               row_colors=cbmicolors_bycbmi,
               col_colors=batchcolors_bycbmi,
               xticklabels=False,
               yticklabels=False,
               figsize=(6,6)
           )
ax = big_ax.ax_heatmap
big_ax.cax.set_visible(False)
N = n_geno.shape[(0)]

ax.axhline(N-39,color="black")
ax.axvline(39,color="black")
ax.text(-4.5,49/2,"OBESE",rotation=90,fontsize=13.5,color="black",verticalalignment="center")
ax.text(-4.5,N-39/2,"LEAN",rotation=90,fontsize=13.5,color="black",verticalalignment="center")
#lt.savefig("./output/replication/our/figures/Figure1c.pdf")
plt.show()


# ### Save batches 1-4 merged
# Finally, we save the merged batches:

n_geno.to_pickle("./output/results/batch1234_geno.p")
pheno.to_pickle("./output/results/batch1234_pheno.p")


# We have thus removed most of the batch effects, and so batches 1 to 4 have been properly merged.
#
# We can now proceed to identify the relevant genes in [this second notebook](2. Extracting the Obesity Signature.ipynb)
