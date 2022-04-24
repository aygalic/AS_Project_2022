library(ggplot2)
library(plotly)  # interactive plots 
library(mvtnorm)
library(rgl)
library(car)

setwd("~/OneDrive/polimi/COURSES/S8/APPLIED_STATS/AS_Project_2022")

#import data
data_patient_= read.delim(file.path("Dataset", "data_clinical_patient.txt"), header = TRUE, comment.char = '#')
data_sample_= read.delim(file.path("Dataset",'data_clinical_sample.txt'), header = TRUE, comment.char = '#')
# use processed data
original_data_mrna_ = read.delim(file.path("Dataset", "1_rpkm.txt"), header = TRUE, comment.char = '#', nrows=5000)



# produce a matrix containing only the genes and cell lines associated with
# a specified cancer type
Build_matrix_for_cancer_type <- function(cancer_name,
                                         data_patient = data_patient_,
                                         data_sample = data_sample_,
                                         original_data_mrna = original_data_mrna_){
  
  #select cells from breast carcinoma patients
  selected_values = data_sample[data_sample$"CANCER_TYPE_DETAILED"==cancer_name,]
  selected_cells = selected_values$SAMPLE_ID
  selected_cells = na.omit(selected_cells)
  indexes_mrna = match(selected_cells, colnames(original_data_mrna))
  mrna_data = original_data_mrna[, c(1, na.omit(indexes_mrna)) ]
  
  #select values only
  data <- as.matrix(mrna_data[,-1])
  
  Y = mrna_data[,1]
  X = colnames(data)
  rownames(data) <- Y
  colnames(data) <- X  
  
  return (data)
}

#   Reorder matrix by how much a gene is expressed
reorder <- function(mat){
  sum = as.data.frame(rowSums(mat))
  return ( as.matrix(scale(mat[order(sum),])) )
}

# Get a table of all cancer types and occurrences
cancer_types <- as.data.frame(table(data_sample_$"CANCER_TYPE_DETAILED"), stringsAsFactors = FALSE)
cancer_types <- cancer_types[order(cancer_types$Freq, decreasing = TRUE),]
names(cancer_types)<-c("Factor", "Freq")


# build matrix for selected cancer types
# The selection of cancer type is done by INDEXES and not by names (1,4,10....)
# Returns a list of the actual matrix and cancer name associated with each observation
#   The matrix is built ordered by cancer type but this is useful to know exactly
#   which cell lines can be grouped together
Build_matrix_for_multiple_cancer_types <- function(selection, types = cancer_types){
  M_ <- NULL
  tag = c()
  
  for(cancer_type in types[selection,1]){
    print(cancer_type)
    m <- Build_matrix_for_cancer_type(cancer_type)
    M_ <- merge(M_, m, all = TRUE, sort = FALSE, by = "row.names")
    row.names(M_) <- M_[,1]
    M_ <- M_[,-1]
    tag <- c(tag, rep(cancer_type, dim(m)[2]))
  }
  
  return(list("Mat" = M_, "Tags" = tag))
}


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

source("src/1_aygalic_pca_generator.R")
reduced_M_scaled = create_reduced_mat(M_scaled, 2)

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



