# This file is meant to handle all the useful function and data structure useful
# for further analysis
#
#

library(ggplot2)
library(plotly)  # interactive plots 
library(mvtnorm)
library(rgl)
library(car)
library(factoextra)

library(htmlwidgets) # saving plot_ly graph
library("tidyverse")





setwd("~/OneDrive/polimi/COURSES/S8/APPLIED_STATS/AS_Project_2022")

#import data
data_patient_= read.delim(file.path("Dataset", "data_clinical_patient.txt"), header = TRUE, comment.char = '#')
data_sample_= read.delim(file.path("Dataset",'data_clinical_sample.txt'), header = TRUE, comment.char = '#')
original_data_mrna_ = read.delim(file.path("Dataset", "1_rpkm.txt"), header = TRUE, comment.char = '#')
data_treatment_auc_ = read.delim(file.path("Dataset", 'data_drug_treatment_auc.txt'), header = TRUE, comment.char = '#')

# preliminary usefull preprocessing
#
# Get a table of all cancer types and occurrences
# this enable us tu have default values available when calling the 
# following functions
cancer_types <- as.data.frame(table(data_sample_$"CANCER_TYPE_DETAILED"), stringsAsFactors = FALSE)
cancer_types <- cancer_types[order(cancer_types$Freq, decreasing = TRUE),]
names(cancer_types)<-c("Factor", "Freq")



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




# build matrix for selected cancer types
# The selection of cancer type is done by INDEXES and not by names (1,4,10....)
# Returns a list of the actual matrix and cancer name associated with each observation
#   The matrix is built ordered by cancer type but this is useful to know exactly
#   which cell lines can be grouped together
Build_matrix_for_multiple_cancer_types <- function(selection = c(1:(dim(cancer_types)-1)),
                                                   types = cancer_types){
  M_ <- NULL
  tag = c()
  
  for(cancer_type in types[selection,1]){
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



## Creating the AUC matrix : code stolen from Robi
create_AUC_matrix <- function(data_patient = data_patient_,
                              data_sample = data_sample_,
                              data_treatment_auc = data_treatment_auc_){
  sample_ID = data_sample$SAMPLE_ID
  
  indexes = match(sample_ID, colnames(data_treatment_auc))

  cancer_data_treatment_auc = data_treatment_auc[, c(na.omit(indexes)) ]
  rownames(cancer_data_treatment_auc) = data_treatment_auc$ENTITY_STABLE_ID
  


  #Remove NA
  data_treatment_auc = cancer_data_treatment_auc
  col_names = data_treatment_auc[1,]
  perc_NA_col = colMeans(is.na(data_treatment_auc))*100
  
  
  row_names = data_treatment_auc[,1]
  perc_NA_row = rowMeans(is.na(data_treatment_auc))*100
  
  max_perc_row = max(perc_NA_row)
  while(max_perc_row > 50)
  {
    index = match(max_perc_row, perc_NA_row)
    data_treatment_auc = data_treatment_auc[-c(index),]
    perc_NA_row = perc_NA_row[-c(index)]
    max_perc_row = max(perc_NA_row)
  }
  
  # 266 --> 218
  
  perc_NA_col = colMeans(is.na(data_treatment_auc))*100
  
  #50% tolerance on the columns
  max_perc_col = max(perc_NA_col)
  while(max_perc_col > 50)
  {
    index = match(max_perc_col, perc_NA_col)
    data_treatment_auc = data_treatment_auc[,-c(index)]
    perc_NA_col = perc_NA_col[-c(index)]
    max_perc_col = max(perc_NA_col)
  }
  
  #from 1065 --> 978
  
  #update the percentages for rows
  perc_NA_row = rowMeans(is.na(data_treatment_auc))*100
  

  max_perc_row = max(perc_NA_row)
  while(max_perc_row > 8)
  {
    index = match(max_perc_row, perc_NA_row)
    data_treatment_auc = data_treatment_auc[-c(index),]
    perc_NA_row = perc_NA_row[-c(index)]
    max_perc_row = max(perc_NA_row)
  }
  
  #dim(data_treatment_auc) #137 rows, 976 columns
  
  #set the mean as value for NA
  rows = rowMeans(data_treatment_auc, na.rm = TRUE)
  
  for(i in 1:137)
  {
    for(j in 1:976)
    {
      if(is.na(data_treatment_auc[i,j])==TRUE)
      {
        data_treatment_auc[i,j] = rows[i]
      }
    }
  }

  cancer_data_treatment_auc = data_treatment_auc
  # Get a table of all cancer types and occurrences
  cancer_types <- as.data.frame(table(data_sample$"CANCER_TYPE_DETAILED"), stringsAsFactors = FALSE)
  cancer_types <- cancer_types[order(cancer_types$Freq, decreasing = TRUE),]
  names(cancer_types)<-c("Factor", "Freq")
  
  cancer_data_treatment_auc = t(cancer_data_treatment_auc)
  
  #I want a new data.frame like the one with AUC value but ordered from the most frequent cancer for better visualization
  
  patients = row.names(cancer_data_treatment_auc)
  M = NULL
  rows = NULL
  for(cancer_type in cancer_types[,1])
  {
    for(p in patients)
    {
      if(cancer_type == data_sample[data_sample$SAMPLE_ID==p, "CANCER_TYPE_DETAILED"])
      {
        M = rbind(M, cancer_data_treatment_auc[p,])
        rows = c(rows, p)
      }
    }
  }
  M = as.data.frame(M)
  rownames(M) = rows
  
  return(M)
}





# Implementation of PCA
create_reduced_mat <- function(mat, n_dim = 2){
  
  # covariance matrix
  cov = cov(mat)
  
  # eigen value decomposition
  decomp = eigen(cov)

  # matrix with the first n_dim eigen vect
  P_reduced = as.matrix(decomp$vectors[,1:n_dim])
  
  # Project the space
  reduced_mat = data.frame(mat%*%P_reduced)
  
  # Make it coherent
  colnames(reduced_mat) <- paste("v", seq(1,n_dim), sep="")
  rownames(reduced_mat) <- rownames(mat)
  
  return(reduced_mat)
}











# We build M with only selected cancer types
# Either provide a matrix of 2 eigen vector
# or the transposed matrix tmat
project_cell_lines_for_cancer_type <- function(selection = c(1:dim(cancer_types)[1]),
                                               space = NULL, tmat = NULL){
  if(is.null(space)){
    if(is.null(mat)){
      print("error, provide either o projection space or a data matrix")
      return()
    }
    
    # PCA as usual
    cov2 = cov(tmat)
    decomp2 = eigen(cov2)
    space = as.matrix(decomp2$vectors[,1:2])
  }
  
  M <- Build_matrix_for_multiple_cancer_types(selection)
  M_ <- M$Mat
  tag <- M$Tags
  
  #scale it
  M_ <- as.matrix(scale(M_))
  tM_ = t(M_)
  
  # We project the cherry picked data on the original space
  reduced_M_ = tM_%*%space
  colnames(reduced_M_) <- c("v1","v2")
  
  plot_ly(data = data.frame(reduced_M_), x = ~v1, y = ~v2,
                 type = "scatter", color = tag) %>%
  layout(
      xaxis = list(range=c(-7,7)),
      yaxis = list(range=c(-9,5)))
}



###
###
###
### Stuff related to the clustering pipeline


contingency_table <- function(k = 2, fun = kmeans, 
                              # param for htrees KEEP THEM NULL IF YOU USE KMEANS
                              hc_method = NULL, hc_metric = NULL) 
{
  
  result.AUC <- NULL
  result.rpkm <- NULL
  
  # case kmeans
  if(is.null(hc_method)){
    result.AUC <- fun(M1, centers=k)$cluster[indexes1]
    result.rpkm <- fun(M2, centers=k)$cluster[indexes2]
  }
  else{
    result.AUC <- fun(M1, k, hc_method = hc_method, hc_metric= hc_metric)$cluster[indexes1]
    result.rpkm <- fun(M2, k, hc_method = hc_method, hc_metric= hc_metric)$cluster[indexes2]
  }
  result <- table(result.AUC, result.rpkm)
  
  #WHAT IS THE AMOUNT OF CORRECT LABELING ?
  result.max <- apply(result, 1, max)
  result.sum <- apply(result, 1, sum)
  
  
  perf <- 100*result.max/result.sum
  
  return(list(result = result, perf = perf, avg_perf=mean(perf)))
}


test_all_cluster_algo <- function(j=10){
  avg_perfs_kmeans <- c()
  
  avg_perfs_ce <- c()
  avg_perfs_ae <- c()
  avg_perfs_se <- c()
  
  avg_perfs_cm <- c()
  avg_perfs_am <- c()
  avg_perfs_sm <- c()
  
  avg_perfs_cc <- c()
  avg_perfs_ac <- c()
  avg_perfs_sc <- c()
  
  for(i in 1:j){avg_perfs_kmeans <- c(avg_perfs_kmeans, contingency_table(i, kmeans)$avg_perf)}
  
  for(i in 1:j){avg_perfs_ce <- c(avg_perfs_ce, contingency_table(i, hcut, hc_method = "complete", hc_metric ="euclidian")$avg_perf)}
  for(i in 1:j){avg_perfs_ae <- c(avg_perfs_ae, contingency_table(i, hcut, hc_method = "average", hc_metric ="euclidian")$avg_perf)}
  for(i in 1:j){avg_perfs_se <- c(avg_perfs_se, contingency_table(i, hcut, hc_method = "single", hc_metric ="euclidian")$avg_perf)}
  
  for(i in 1:j){avg_perfs_cm <- c(avg_perfs_cm, contingency_table(i, hcut, hc_method = "complete", hc_metric ="manhattan")$avg_perf)}
  for(i in 1:j){avg_perfs_am <- c(avg_perfs_am, contingency_table(i, hcut, hc_method = "average", hc_metric ="manhattan")$avg_perf)}
  for(i in 1:j){avg_perfs_sm <- c(avg_perfs_sm, contingency_table(i, hcut, hc_method = "single", hc_metric ="manhattan")$avg_perf)}
  
  for(i in 1:j){avg_perfs_cc <- c(avg_perfs_cc, contingency_table(i, hcut, hc_method = "complete", hc_metric ="canberra")$avg_perf)}
  for(i in 1:j){avg_perfs_ac <- c(avg_perfs_ac, contingency_table(i, hcut, hc_method = "average", hc_metric ="canberra")$avg_perf)}
  for(i in 1:j){avg_perfs_sc <- c(avg_perfs_sc, contingency_table(i, hcut, hc_method = "single", hc_metric ="canberra")$avg_perf)}
  
  
  
  
  n_clusters <- c(1:j)
  
  result <- data.frame(n_clusters, avg_perfs_kmeans,
                       avg_perfs_ce, avg_perfs_ae, avg_perfs_se,
                       avg_perfs_cm, avg_perfs_am, avg_perfs_sm,
                       avg_perfs_cc, avg_perfs_ac, avg_perfs_sc)
  
  
  fig <- plot_ly(result, x = ~n_clusters, y = ~avg_perfs_kmeans, name = 'kmeans', type = 'scatter', mode = 'lines+markers') 
  fig <- fig %>% add_trace(y = ~avg_perfs_ce, name = 'complete euclidian', mode = 'lines+markers') 
  fig <- fig %>% add_trace(y = ~avg_perfs_ae, name = 'average euclidian', mode = 'lines+markers') 
  fig <- fig %>% add_trace(y = ~avg_perfs_se, name = 'single euclidian', mode = 'lines+markers') 
  
  fig <- fig %>% add_trace(y = ~avg_perfs_cm, name = 'complete manhattan', mode = 'lines+markers') 
  fig <- fig %>% add_trace(y = ~avg_perfs_am, name = 'average manhattan', mode = 'lines+markers') 
  fig <- fig %>% add_trace(y = ~avg_perfs_sm, name = 'single manhattan', mode = 'lines+markers') 
  
  fig <- fig %>% add_trace(y = ~avg_perfs_cc, name = 'complete canberra"', mode = 'lines+markers') 
  fig <- fig %>% add_trace(y = ~avg_perfs_ac, name = 'average canberra', mode = 'lines+markers') 
  fig <- fig %>% add_trace(y = ~avg_perfs_sc, name = 'single canberra', mode = 'lines+markers') 
  
  return(list(result = result, fig = fig))
} 


# This function is designed to match the centroid of the 2 groups of clusters
# in a given projection space
#
# Match clusters2 to clusters1
#
# returns a new list of clusters based on clusters2, 
# but with indexes swapped to match clusters1's centroids
align_clusters_in_space <- function(clusters1, clusters2, projected_obs){
  n_clusters <- max(clusters1)
  if(n_clusters<2) return(clusters2)
  
  new_clusters <- clusters2
  for(i in 1:n_clusters){
    # find centroid of each cluster for the first dataset
    data <- projected_obs[clusters1==i,]
    centroid1 <- rapply(data, mean)
    
    # find the matching centroid
    # start with clust 1 ofc
    k = 1
    centroid2 <- rapply(projected_obs[clusters2==k,], mean)
    # compute the original distance
    d = dist(rbind(centroid1, centroid2))
    
    
    for(j in 2:n_clusters){
      centroid2_ <- rapply(projected_obs[clusters2==j,], mean)
      # compute the distance bewteen current centroids
      d_ = dist(rbind(centroid1, centroid2_))
      if(d_ < d){
        d <- d_
        centroid2 <- centroid2_
        k = j
      }
    }
    # set the new cluster k of clusters2 to the cluster i of clusters1
    new_clusters[clusters2==k] <- i
  }
  return (new_clusters)
}


