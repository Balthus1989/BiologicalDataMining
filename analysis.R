library(ggfortify)
library(limma)

setwd("/home/gianmarco//Scrivania/First_Semester/Data_Mining_Laboratory/p/DataMiningProject")

batchTot <- read.table("./output/Ztot.txt", sep='\t', header=TRUE, row.names = 1)

batchTot[,'Cbmi'] <- as.factor(batchTot[,'Cbmi'])
batchTot[,'batches'] <- as.factor(batchTot[,'batches'])

Z.pca <- prcomp(batchTot, cente=TRUE, scale.=TRUE)
print(Z.pca)
plot(Z.pca, type='l')
summary(Z.pca)
autoplot(Z.pca, data = batchTot, colour = 'Cbmi')
autoplot(Z.pca, data = batchTot, colour = 'batches')

Z.1 <- removeBatchEffect(batchTot[,!(names(batchTot)%in% c("Cbmi", "batches"))])
Z1.pca <- prcomp(Z.1, cente=TRUE, scale.=TRUE)
print(Z1.pca)
plot(Z1.pca, type='l')
summary(Z1.pca)
autoplot(Z1.pca, data = batchTot, colour = 'Cbmi')
autoplot(Z1.pca, data = batchTot, colour = 'batches')