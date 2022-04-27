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
original_data_mrna_ = read.delim(file.path("Dataset", "1_rpkm.txt"), header = TRUE, comment.char = '#')

# Get a table of all cancer types and occurrences
cancer_types <- as.data.frame(table(data_sample_$"CANCER_TYPE_DETAILED"), stringsAsFactors = FALSE)
cancer_types <- cancer_types[order(cancer_types$Freq, decreasing = TRUE),]
names(cancer_types)<-c("Factor", "Freq")


# creating a matrix with all cancer types
#
# We could also create a matrix with specified cancer types only
cancer_type_selection <- c(1:(length(cancer_types$Factor)-1))
M = Build_matrix_for_multiple_cancer_types(cancer_type_selection)
tags = M$Tags
M <- M$Mat
M_scaled <- as.matrix(scale(M))

M <- t(M_scaled)

# k-means
result.k <- kmeans(M, centers=3) # Centers: fixed number of clusters




# WHAT HAPPENS IF YOU PLOT THOSES ON THE REDUCED SPACE PROVIDED BY PCA ?
# create a projected space 
reduced_M_scaled = create_reduced_mat(M, 2)


fviz_nbclust(M, FUN = kmeans, method = "silhouette") 
fviz_nbclust(M, FUN = kmeans, method = "wss")
# 3 clusters being optimal, could go to 7



plot_ly(data = data.frame(reduced_M_scaled), x = ~v1, y = ~v2,
        text = rownames(reduced_M_scaled), type = "scatter", 
        color = factor(result.k$cluster)) %>% 
  layout(margin = c(10,10,10,10,0))





# experimenting with the number of clusters
result.k <- kmeans(M, centers=7) # Centers: fixed number of clusters
plot_ly(data = data.frame(reduced_M_scaled), x = ~v1, y = ~v2,
        text = rownames(reduced_M_scaled), type = "scatter", 
        color = factor(result.k$cluster)) %>% 
  layout(margin = c(10,10,10,10,0))

# experimenting with the number of clusters
result.k <- kmeans(M, centers=5) # Centers: fixed number of clusters
plot_ly(data = data.frame(reduced_M_scaled), x = ~v1, y = ~v2,
        text = rownames(reduced_M_scaled), type = "scatter", 
        color = factor(result.k$cluster)) %>% 
  layout(margin = c(10,10,10,10,0))



