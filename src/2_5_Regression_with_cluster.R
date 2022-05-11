setwd("~/OneDrive/polimi/COURSES/S8/APPLIED_STATS/AS_Project_2022")

source("src/utilities.R")

################################################
############### PREPARE THE DATA ###############
################################################

# steal Robi's work
M_AUC <- create_AUC_matrix()
M1<- as.matrix(scale(M_AUC))

# merge with mine
M <- Build_matrix_for_multiple_cancer_types()$Mat
M_scaled <- as.matrix(scale(M))
M2 <- M_scaled


# select row from M1 that are present in M2
indexes1 <- rownames(M1) %in% colnames(M2)
M1_ <- M1 [indexes1,]
dim(M1_)

# select cols from M2 that are present in M1
indexes2 <- colnames(M2) %in% rownames(M1)
M2_ <- M2[,indexes2]
dim(M2_)

# the selection is a little different than before since are going to transpose 
# the matrixlater




#################################################
################# ORIGINAL TEST #################
#################################################

# Do k_means on original matrix and keep only the labels of interest
k=10
result.rpkm.genes <- kmeans(M2, centers=k)$cluster


# NOW WE REDUCE THE CELL LINE TO THEIR RESPECTIVE CLUSTERS
result.rpkm.genes



# Can we visualize what we are doing with PCA ?
M2_reduced <- create_reduced_mat(M2_)

plot_ly(data = data.frame(M2_reduced), x = ~v1, y = ~v2,
        text = rownames(M2_reduced), type = "scatter",
        color = factor(result.rpkm.genes)) %>% 
  layout(margin = c(10,10,10,10,0))


# NOW : reduce each gene to the average of all it's clusters :

plot_ly(data = data.frame(M2_reduced), x = ~v1, y = ~v2,
        text = rownames(M2_reduced), type = "scatter",
        color = factor(result.rpkm.genes)) %>% 
  layout(margin = c(10,10,10,10,0))



plot_ly(data = data.frame(M2_reduced[result.rpkm.genes==1,]), x = ~v1, y = ~v2,
        text = rownames(M2_reduced[result.rpkm.genes==1,]), type = "scatter") %>% 
  layout(margin = c(10,10,10,10,0))


# We first select genes of cluster 1, make the mean
M_compressed <- data.frame("cluster" = apply(as.data.frame(M2_[result.rpkm.genes==1,]), 2, mean))
# then we do it for all remaning clusters
for(i in 2:k){
  M_compressed <- data.frame(M_compressed, "cluster" = apply(as.data.frame(M2_[result.rpkm.genes==i,]), 2, mean))
}


# try to do regression

X <- M_compressed
X$y <- M1_[,1]

fit <- lm(y ~ ., X)
summary(fit)




