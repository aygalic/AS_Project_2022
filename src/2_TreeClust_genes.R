library(ggplot2)
library(plotly)  # interactive plots 
library(mvtnorm)
library(rgl)
library(car)

setwd("~/OneDrive/polimi/COURSES/S8/APPLIED_STATS/AS_Project_2022")
# all the useful functions for pca and projections
source("src/utilities.R")



#import data
data_patient_= read.delim(file.path("Dataset", "data_clinical_patient.txt"), header = TRUE, comment.char = '#')
data_sample_= read.delim(file.path("Dataset",'data_clinical_sample.txt'), header = TRUE, comment.char = '#')
# use processed data
# original_data_mrna_ = read.delim(file.path("Dataset", "1_rpkm.txt"), header = TRUE, comment.char = '#', nrows=5000)
# we might as well use all the data
original_data_mrna_ = read.delim(file.path("Dataset", "1_rpkm.txt"), header = TRUE, comment.char = '#')





# Get a table of all cancer types and occurrences
cancer_types <- as.data.frame(table(data_sample_$"CANCER_TYPE_DETAILED"), stringsAsFactors = FALSE)
cancer_types <- cancer_types[order(cancer_types$Freq, decreasing = TRUE),]
names(cancer_types)<-c("Factor", "Freq")




# creating a matrix will all cancer types
#
# We could also create a matrix with specified cancer types only
cancer_type_selection <- c(1:(length(cancer_types$Factor)-1))
M = Build_matrix_for_multiple_cancer_types(cancer_type_selection)
tags = M$Tags
M <- M$Mat
M_scaled <- as.matrix(scale(M))



M_ <- M_scaled

M_.e <- dist(M_, method='euclidean')
M_.m <- dist(M_, method='manhattan')
M_.c <- dist(M_, method='canberra')




# making every possible tree
M_.es <- hclust(M_.e, method='single')
M_.ea <- hclust(M_.e, method='average')
M_.ec <- hclust(M_.e, method='complete')

M_.ms <- hclust(M_.m, method='single')
M_.ma <- hclust(M_.m, method='average')
M_.mc <- hclust(M_.m, method='complete')

M_.cs <- hclust(M_.c, method='single')
M_.ca <- hclust(M_.c, method='average')
M_.cc <- hclust(M_.c, method='complete')



# if we want more detailed information on euclidean-complete
# clustering:
names(M_.ec)
M_.ec$merge  # order of aggregation of statistical units / clusters
M_.ec$height # distance at which we have aggregations
M_.ec$order  # ordering that allows to avoid intersections in the dendrogram

# plot of the dendrograms
x11()
par(mfrow=c(1,3))
plot(M_.es, main='euclidean-single', hang=-0.1, xlab='', labels=F, cex=0.6, sub='')
plot(M_.ec, main='euclidean-complete', hang=-0.1, xlab='', labels=F, cex=0.6, sub='')
plot(M_.ea, main='euclidean-average', hang=-0.1, xlab='', labels=F, cex=0.6, sub='')

dev.off()





# plot dendrograms (2 clusters)
x11()
par(mfrow=c(1,3))
plot(M_.es, main='euclidean-single', hang=-0.1, xlab='', labels=F, cex=0.6, sub='')
rect.hclust(M_.es, k=2)
plot(M_.ec, main='euclidean-complete', hang=-0.1, xlab='', labels=F, cex=0.6, sub='')
rect.hclust(M_.ec, k=2)
plot(M_.ea, main='euclidean-average', hang=-0.1, xlab='', labels=F, cex=0.6, sub='')
rect.hclust(M_.ea, k=2)



# How to cut a dendrogram?
# We generate vectors of labels through the command cutree()
help(cutree)

# Fix k=2 clusters:
cluster.ec <- cutree(M_.ec, k=2) # euclidean-complete:
cluster.ec




# WHAT HAPPENS IF YOU PLOT THOSES ON THE REDUCED SPACE PROVIDED BY PCA ?
# create a projected space 
reduced_M_scaled = create_reduced_mat(M_scaled, 2)

# big brain graph
k=3
cluster.es <- cutree(M_.es, k=k)
cluster.ea <- cutree(M_.ea, k=k)
cluster.ec <- cutree(M_.ec, k=k)

cluster.ms <- cutree(M_.ms, k=k)
cluster.ma <- cutree(M_.ma, k=k)
cluster.mc <- cutree(M_.mc, k=k)

cluster.cs <- cutree(M_.cs, k=k)
cluster.ca <- cutree(M_.ca, k=k)
cluster.cc <- cutree(M_.cc, k=k)

x11()
par(mfrow=c(3,3))
plot(data.frame(reduced_M_scaled), main = 'euclidean single',   col=cluster.es, pch=16, asp=1)
plot(data.frame(reduced_M_scaled), main = 'euclidean average',  col=cluster.ea, pch=16, asp=1)
plot(data.frame(reduced_M_scaled), main = 'euclidean complete', col=cluster.ec, pch=16, asp=1)
plot(data.frame(reduced_M_scaled), main = 'manhattan single',   col=cluster.ms, pch=16, asp=1)
plot(data.frame(reduced_M_scaled), main = 'manhattan average',  col=cluster.ma, pch=16, asp=1)
plot(data.frame(reduced_M_scaled), main = 'manhattan complete', col=cluster.mc, pch=16, asp=1)
plot(data.frame(reduced_M_scaled), main = 'canberra single',    col=cluster.cs, pch=16, asp=1)
plot(data.frame(reduced_M_scaled), main = 'canberra average',   col=cluster.ca, pch=16, asp=1)
plot(data.frame(reduced_M_scaled), main = 'canberra complete',  col=cluster.cc, pch=16, asp=1)
dev.off()


