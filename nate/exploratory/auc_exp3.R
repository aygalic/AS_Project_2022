#attempt at clustering drugs for a specific cancer type

source("../utils/nate_utils.R")
source("../utils/exp_base_script.R")

cancer.data<-block_dat("CENTRAL_NERVOUS_SYSTEM", auc)


#lets see if there's some non sparse cell lines 
cells.sparcity = sort(colMeans(is.na(cancer.data)), decreasing = FALSE)
cells.sparcity

#note: the following choice might not reflect genetic diversity in patients, chosen
#simply since it would result in least amount of information lost 
cells.viable = c(names(cells.sparcity)[1:5])
viable.data = cancer.data[,cells.viable]


#remove some NAs (this also removes some drugs though!)
dim(viable.data)
viable.data = na.omit(viable.data)
dim(viable.data)


#quick plot to see intuitively if some groups 
viable.avg = rowMeans(viable.data)
x11()
plot(viable.avg, rep(0, length(viable.avg)), pch = 19, xlab="avg auc")

#cluster into what we suspect is 2 clusters for drugs 
#clustering done over p dimensions instead of average (robustness of model)
drugs.dist = dist(viable.data, method = "euclidean")
drugs.hclust = hclust(drugs.dist, method = "complete")

drugs.clust = cutree(drugs.hclust, k=2)

#dendrogram 
x11()
plot(drugs.hclust)

#pairwise plots
x11()
plot(viable.data, col=drugs.clust+1, pch = 19)

#2D plot 
x11()
plot(viable.avg, rep(0, length(viable.avg)), col = drugs.clust+1, 
     pch = 19, xlab = "Average AUC", ylab = "", main="AUC of drugs on cancer type")


#test for gaussianity between clusters of drugs 
load("/Users/nathanielnethercott/Desktop/School/Polimi/2021:2022/AS/Labs-20220417/Lab 5/mcshapiro.test.RData")

#In both cases we don't have gaussianity ://
#original variables
mcshapiro.test(viable.data[drugs.clust == '1',])$pvalue
mcshapiro.test(viable.data[drugs.clust == '2',])$pvalue

#averaged 
g1.pval = shapiro.test(viable.avg[drugs.clust == '1'])$p
g2.pval = shapiro.test(viable.avg[drugs.clust == '2'])$p

x11()
par(mfrow=(c(1,2)))
hist(viable.avg[drugs.clust == '1'], main = paste("hist for group 1, p=", round(g1.pval,4)), xlim = (0:1), xlab = "avg auc")
hist(viable.avg[drugs.clust == '2'], main = paste("hist for group 2, p=", round(g2.pval,4)), xlim = (0:1), xlab = "avg auc")


#TODO:
  #- reference the drugs_labelled dataset to see if any underlying traits between clusters 
  #- manually check to make sure averaging and stuff done right 

