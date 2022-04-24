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
original_data_mrna_ = read.delim(file.path("Dataset", "1_rpkm.txt"), header = TRUE, comment.char = '#', nrows=5000)





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


M.e <- dist(M, method='euclidean')
M.m <- dist(M, method='manhattan')
M.c <- dist(M, method='canberra')


# actually, the data are never ordered according to (unknown) labels
n=100
misc <- sample(n)
M_ <- M[misc,]

# Bypass the sampling
M_ <- M_scaled

M_.e <- dist(M_, method='euclidean')
M_.m <- dist(M_, method='manhattan')
M_.c <- dist(M_, method='canberra')

x11()
image(1:n,1:n,as.matrix(M_.e), main='metrics: Euclidean', asp=1, xlab='i', ylab='j' )
graphics.off()



M_.es <- hclust(M_.e, method='single')
M_.ea <- hclust(M_.e, method='average')
M_.ec <- hclust(M_.e, method='complete')



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

cluster.ec <- cutree(M_.ec, k=3)

plot_ly(data = data.frame(reduced_M_scaled), x = ~v1, y = ~v2,
        text = rownames(reduced_M_scaled), type = "scatter", 
        color = cluster.ec) %>% 
  layout(margin = c(10,10,10,10,0))





# experimenting with the number of clusters
cluster.ec <- cutree(M_.ec, k=3)
plot_ly(data = data.frame(reduced_M_scaled), x = ~v1, y = ~v2,
        text = rownames(reduced_M_scaled), type = "scatter", 
        color = factor(cluster.ec)) %>% 
  layout(margin = c(10,10,10,10,0))

# experimenting with the number of clusters
cluster.ec <- cutree(M_.ec, k=10)
plot_ly(data = data.frame(reduced_M_scaled), x = ~v1, y = ~v2,
        text = rownames(reduced_M_scaled), type = "scatter", 
        color = factor(cluster.ec)) %>% 
  layout(margin = c(10,10,10,10,0))



