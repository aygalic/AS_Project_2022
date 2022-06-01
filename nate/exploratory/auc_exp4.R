source("../utils/nate_utils.R")
source("../utils/exp_base_script.R")

#let's do clustering on two cancers w two groups to see something
#eventually scale this process to a plot with multiple cancer pairings 

sort(table(cancer_freqs), decreasing = TRUE)
cancers = c("BREAST", "CENTRAL_NERVOUS_SYSTEM")

cancer.data = block_dat(cancers, auc)

library(zoo)
cancer.data = na.aggregate(cancer.data)

#get samples as auc scores 
cancer.data = t(cancer.data)

### PCA
cancer.cov = cov(cancer.data)
cancer.eigen = eigen(cancer.cov)
scores = cancer.data%*%cancer.eigen$vectors

#clustering on un-projected
cancer.dist = dist(cancer.data, method="euclidean")
cancer.hclust = hclust(cancer.dist, method = "complete")

x11()
plot(cancer.hclust)

cluster.ec <- cutree(cancer.hclust, k=2) 

#now get `real` labels and colour maps
labels.real = str_util(rownames(cancer.data))
reference_map = ifelse(labels.real == "BREAST", 'red', 'blue')
cluster_map = ifelse(cluster.ec==1,'red','blue')

x11()
par(mfrow=c(1,2))
plot(scores[,1], scores[,2], col=cluster_map, pch=16, asp=1, main="hierarchical clustering", lwd=2)
plot(scores[,1], scores[,2], col=reference_map, pch=16, asp=1, main="real labels", lwd=2)


#quantify the results
metric = data.frame(cancer = labels.real, label = cluster.ec)
table(metric)


#clustering on pc's

#let's investigate the principle components some more 
x11()
par(mfrow=c(1,2))
plot(cancer.eigen$values, xlab="pc", ylab="variance")
plot(cumsum(cancer.eigen$values)/sum(cancer.eigen$values), xlab="pc", ylab="explained variance")
abline(h=0.85, col='red', lty=2)

x11()
par(mfrow=c(3,1))
barplot(cancer.eigen$vectors[,1])
barplot(cancer.eigen$vectors[,2])
barplot(cancer.eigen$vectors[,3])

#select number of principle components corresponding to 90% of explained variance 
n.comp = min(which(cumsum(cancer.eigen$values)/sum(cancer.eigen$values)>=0.9))
n.comp

#now cluster again 
cancer.reduced = scores[,1:n.comp]
reduced.dist = dist(cancer.reduced, method="euclidean")
reduced.hclust = hclust(reduced.dist, method = "complete")

x11()
plot(reduced.hclust)

reduced.ec <- cutree(reduced.hclust, k=2) 

#plot
cluster_map.reduced = ifelse(reduced.ec==1,'red','blue')

x11()
par(mfrow=c(1,2))
plot(scores[,1], scores[,2], col=cluster_map.reduced, pch=16, asp=1, main="hierarchical clustering", lwd=2)
plot(scores[,1], scores[,2], col=reference_map, pch=16, asp=1, main="real labels", lwd=2)

metric.reduced = data.frame(cancer = labels.real, label = reduced.ec)
table(metric.reduced)

#-------------------------------------------------------------------------------
#for luca's app 
#tosave = data.frame(cancer.data)
#save(tosave, file = "breast_cns.Rda")

#-------------------------------------------------------------------------------
#Now let's do the same type of thing with fewer treatments 
cancer2.data = block_dat(cancers, auc)
na_rows = sort(rowSums(is.na(cancer2.data)), decreasing=TRUE)

#we can choose to reduce num rows such that ratio of nas to sample size is larger than 50%
num_samples = dim(cancer2.data)[2]
drugs_kept = names(which(na_rows/num_samples <0.5))

#reduce 
cancer2.data = cancer2.data[drugs_kept, ]

#go through the clustering once again 
cancer2.data = na.aggregate(cancer2.data)
cancer2.data = t(cancer2.data)

### PCA
cancer2.cov = cov(cancer2.data)
cancer2.eigen = eigen(cancer2.cov)
scores2 = cancer2.data%*%cancer2.eigen$vectors

#clustering on un-projected
cancer2.dist = dist(cancer2.data, method="euclidean")
cancer2.hclust = hclust(cancer2.dist, method = "complete")

x11()
plot(cancer2.hclust)

cluster2.ec <- cutree(cancer2.hclust, k=2) 

#now get `real` labels and colour maps
labels2.real = str_util(rownames(cancer2.data))
reference_map2 = ifelse(labels2.real == "BREAST", 'red', 'blue')
cluster_map2 = ifelse(cluster2.ec==1,'red','blue')

x11()
par(mfrow=c(1,2))
plot(scores2[,1], scores2[,2], col=cluster_map2, pch=16, asp=1, main="hierarchical clustering", lwd=2)
plot(scores2[,1], scores2[,2], col=reference_map2, pch=16, asp=1, main="real labels", lwd=2)

metric2 = data.frame(cancer = labels2.real, label = cluster2.ec)
table(metric2)

