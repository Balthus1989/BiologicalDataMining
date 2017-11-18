import pandas as pd

their_totData = pd.read_pickle("./obesity-score-master/output/batches1-4_merged/batch1234_geno.p")
#print(their_totData)
our_totData = pd.read_pickle("./output/replication/our/batches1-4_merged/batch1234_geno.p")
#print(our_totData)

their_sigs = pd.read_pickle("./obesity-score-master/output/signature/signature.p")
our_sigs = pd.read_pickle("./output/replication/our/signature/signature.p")

#print(their_sigs)
#print(our_sigs)

with open("./output/replication/final_comparison.txt", 'w') as handle:
    handle.write("IntersectionSigs: "+str(len(set.intersection(set(their_sigs.index), set(our_sigs.index))))\
    +"\nJaccardDist: "+str(len(set.intersection(set(their_sigs.index), set(our_sigs.index)))/len(set.union(set(their_sigs.index), set(our_sigs.index))))
    +"\nIntersection: "+str(len(set.intersection(set(our_totData.columns), set(their_totData.columns))))\
    +"\n#TheirGenes: "+str(len(their_totData.columns))\
    +"\n#OurGenes: "+str(len(our_totData.columns))\
    +"\nJaccardDist: "+str(len(set.intersection(set(their_totData.index), set(our_totData.index)))/len(set.union(set(their_totData.index), set(our_totData.index)))))

print("IntersectionSigs: "+str(len(set.intersection(set(their_sigs.index), set(our_sigs.index))))\
+"JaccardDist: "+str(len(set.intersection(set(their_sigs.index), set(our_sigs.index)))/len(set.union(set(their_sigs.index), set(our_sigs.index)))))
print("\nIntersection: "+str(len(set.intersection(set(our_totData.columns), set(their_totData.columns))))\
+"\n#TheirGenes: "+str(len(their_totData.columns))\
+"\n#OurGenes: "+str(len(our_totData.columns))\
+"\nJaccardDist: "+str(len(set.intersection(set(their_totData.index), set(our_totData.index)))/len(set.union(set(their_totData.index), set(our_totData.index)))))
