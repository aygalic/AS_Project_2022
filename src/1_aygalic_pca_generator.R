library(ggplot2)
library(plotly)  # interactive plots 

setwd("~/OneDrive/polimi/COURSES/S8/APPLIED_STATS/AS_Project_2022")

#import data
data_patient_= read.delim(file.path("Dataset", "data_clinical_patient.txt"), header = TRUE, comment.char = '#')
data_sample_= read.delim(file.path("Dataset",'data_clinical_sample.txt'), header = TRUE, comment.char = '#')
original_data_mrna_ = read.delim(file.path("Dataset", "data_mrna_seq_rpkm.txt"), header = TRUE, comment.char = '#', nrows=5000)
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

# Get a table of all cancer types and occurence
cancer_types <- as.data.frame(table(data_sample_$"CANCER_TYPE_DETAILED"), stringsAsFactors = FALSE)
cancer_types <- cancer_types[order(cancer_types$Freq, decreasing = TRUE),]
names(cancer_types)<-c("Factor", "Freq")



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

###########################################
############ Building Heatmaps ############ 
###########################################

# creating a matrix will all cancer types
#
# We could also create a matrix with specified cancer types only
cancer_type_selection <- c(1:(length(cancer_types$Factor)-1))
M = Build_matrix_for_multiple_cancer_types(cancer_type_selection)
tags = M$Tags
M <- M$Mat
M_scaled <- as.matrix(scale(M))


# plot heatmap
plot_ly(x=colnames(M_scaled), y=rownames(M_scaled), z = M_scaled, type = "heatmap")





#############################
############ PCA ############
#############################

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


# We reduce the space of cell line and project genes (to see if there are any gene of interest ?)

reduced_M_scaled = create_reduced_mat(M_scaled, 2)

plot_ly(data = data.frame(reduced_M_scaled), x = ~v1, y = ~v2,
               text = rownames(reduced_M_scaled), type = "scatter") %>% 
               layout(margin = c(10,10,10,10,0))



# adding color based on overall gene expression (before scaling)

gene_exp = apply(M, 1, sum)

plot_ly(data = data.frame(reduced_M_scaled), x = ~v1, y = ~v2,
        text = rownames(reduced_M_scaled), type = "scatter",
        color = gene_exp) %>% 
  layout(margin = c(10,10,10,10,0))


# same but log scaled ?

gene_exp = apply(M, 1, sum)
gene_exp_scaled = log(gene_exp)

plot_ly(data = data.frame(reduced_M_scaled), x = ~v1, y = ~v2,
        text = rownames(reduced_M_scaled), type = "scatter",
        color = gene_exp_scaled) %>% 
  layout(margin = c(10,10,10,10,0))



# 3D plot 

reduced_M_scaled = create_reduced_mat(M_scaled, 3)

plot_ly(data = data.frame(reduced_M_scaled), x = ~v1, y = ~v2, z = ~v3,
        text = rownames(reduced_M_scaled),
        color = gene_exp_scaled) %>% 
  layout(margin = c(10,10,10,10,0))









# We want to do PCA projecting the cell lines on a reduced gene space
tM_scaled = t(M_scaled)

cov2 = cov(tM_scaled)
decomp2 = eigen(cov2)
barplot(decomp2$values[1:50]/sum(decomp2$values), ylim = c(0,1))

P2 = as.matrix(decomp2$vectors[,1:2])



# We build M with only selected cancer types
project_cell_lines_for_cancer_type <- function(selection = c(1:dim(cancer_types)[1]),
                                               space = P2){
  
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


# next thing to do :
# Do PCA while hiding some of the classes from the data :
# only plot breast & lung or any combinaison

# plot all of them
project_cell_lines_for_cancer_type()

#plot only breast cancer
project_cell_lines_for_cancer_type(c(7))

#plot  breast cancer & a few other
project_cell_lines_for_cancer_type(c(11,19,15,7))
project_cell_lines_for_cancer_type(c(13,15,7))


# quick pca on gene deprived matrix
reduced_M_dep = create_reduced_mat(M_gene_deprived, 2)


plot_ly(data = data.frame(reduced_M_dep), x = ~v1, y = ~v2,
               text = rownames(reduced_M_dep), type = "scatter") %>% 
  layout(margin = c(10,10,10,10,0))




