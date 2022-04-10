library(ggplot2)
library(plotly)  # interactive plots 

setwd("~/OneDrive/polimi/COURSES/S8/APPLIED_STATS/AS_Project_2022")


#import data
data_patient_= read.delim(file.path("Dataset", "data_clinical_patient.txt"), header = TRUE, comment.char = '#')
data_sample_= read.delim(file.path("Dataset",'data_clinical_sample.txt'), header = TRUE, comment.char = '#')
original_data_mrna_ = read.delim(file.path("Dataset", "data_mrna_seq_rpkm.txt"), header = TRUE, comment.char = '#', nrows=100)


Bmatrix_for_cancer_type <- function(cancer_name,
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




cancer_types <- as.data.frame(table(data_sample_$"CANCER_TYPE_DETAILED"), stringsAsFactors = FALSE)
cancer_types <- cancer_types[order(cancer_types$Freq, decreasing = TRUE),]
names(cancer_types)<-c("Factor", "Freq")

m1 <- Bmatrix_for_cancer_type(cancer_types$Factor[1])


# creating a matrix will all cancer types
M = m1
for(cancer_type in cancer_types$Factor[-length(cancer_types$Factor)]){
  print(cancer_type)
  m <- Bmatrix_for_cancer_type(cancer_type)
  M <- merge(M,m, all = TRUE, sort = FALSE, by = "row.names")
  row.names(M) <- M[,1]
  M <- M[,-1]
}


###########################
###########################   Let's maybe not do that
###########################
sum = as.data.frame(rowSums(M))
length(order(sum))
dim(M)
M[order(sum),]
#dev.off()
plot.default(sum[order(sum),])
#order
#M <- as.matrix(scale(M[order(sum),]))
###########################
###########################
###########################



M <- as.matrix(scale(M))





Y = rownames(M)
X = colnames(M)
fig <- plot_ly(x=X, y=Y, z = M, type = "heatmap")
fig

# Dendrogram
library(heatmaply)
heatmaply(M, k_col = length(cancer_types$Factor) -1)


#############################
############ PCA ############
#############################

# We reduce the space of cell line and project genes (to see if there are any gene of interet ?)

# covariance matrix
cov = cov(M)

# eigen value decomposition
decomp = eigen(cov)


# ploting the first loadings
barplot(decomp$values[1:10]/sum(decomp$values), ylim = c(0,1))


# matrix with the first 2 eigen vect
P = as.matrix(decomp$vectors[,1:2])
# not sure how usefull is that
reduced_M = data.frame(M%*%P)



colnames(reduced_M) <- c("v1","v2")
rownames(reduced_M) <- rownames(M)
plot(reduced_M)



fig <- plot_ly(data = data.frame(reduced_M), x = ~v1, y = ~v2,
               text = rownames(reduced_M), type = "scatter") %>% 
               layout(margin = c(10,10,10,10,0))


fig

# Double check that PCA is done correctly




# We want to do PCA projecting the cell lines on a reduced gene space
tM = t(M)



cov2 = cov(tM)
decomp2 = eigen(cov2)
barplot(decomp2$values[1:50]/sum(decomp2$values), ylim = c(0,1))

P2 = as.matrix(decomp2$vectors[,1:2])
barplot(P2[1:50,1])

reduced_M = tM%*%P2



colnames(reduced_M) <- c("v1","v2")
plot(reduced_M)
fig <- plot_ly(data = data.frame(reduced_M), x = ~v1, y = ~v2, type = "scatter")

fig


tag = rownames(reduced_M)
for(i in 1:(length(tag)-1)){
  tag[[i]] <- paste(strsplit(tag[[i]], "_")[[1]][-1],  collapse = " ")
}
tag

fig <- plot_ly(data = data.frame(reduced_M), x = ~v1, y = ~v2, type = "scatter", color = tag)

fig

# next thing to do :
# try to do pca while hiding some of the classes from the data :
# only plot breast & lung or any combinaison










# try to remove the most expressed genes from the heatmap to see if it become readable
#
#