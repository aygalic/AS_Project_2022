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
big_boy_matrix = m1
for(cancer_type in cancer_types$Factor[-length(cancer_types$Factor)]){
  print(cancer_type)
  m <- Bmatrix_for_cancer_type(cancer_type)
  big_boy_matrix <- merge(big_boy_matrix,m, all = TRUE, sort = FALSE, by = "row.names")
  row.names(big_boy_matrix) <- big_boy_matrix[,1]
  big_boy_matrix <- big_boy_matrix[,-1]
}




sum = as.data.frame(rowSums(big_boy_matrix))
length(order(sum))
dim(big_boy_matrix)
big_boy_matrix[order(sum),]
#dev.off()
plot.default(sum[order(sum),])
#order
#big_boy_matrix <- as.matrix(scale(big_boy_matrix[order(sum),]))




big_boy_matrix <- as.matrix(scale(big_boy_matrix))





Y = rownames(big_boy_matrix)
X = colnames(big_boy_matrix)
fig <- plot_ly(x=X, y=Y, z = big_boy_matrix, type = "heatmap")
fig

library(heatmaply)
heatmaply(big_boy_matrix, k_col = length(cancer_types$Factor) -1)


#############################
############ PCA ############
#############################

# We want to do PCA on the cell lines first
cov = cov(big_boy_matrix)

decomp = eigen(cov)

barplot(decomp$values[1:10]/sum(decomp$values), ylim = c(0,1))

# matrix with the first 2 eigen vect
P = as.matrix(decomp$vectors[,1:2])
# not sure how usefull is that
barplot(P[1:50,1])
reduced_P = big_boy_matrix%*%P



colnames(reduced_P) <- c("v1","v2")
plot(reduced_P)
fig <- plot_ly(data = data.frame(reduced_P), x = ~v1, y = ~v2, type = "scatter")

fig


# We want to do PCA on the genes
t = t(big_boy_matrix)



cov2 = cov(t)
decomp2 = eigen(cov2)
barplot(decomp2$values[1:50]/sum(decomp2$values), ylim = c(0,1))

P2 = as.matrix(decomp2$vectors[,1:2])
barplot(P2[1:50,1])

reduced_P = t%*%P2



colnames(reduced_P) <- c("v1","v2")
plot(reduced_P)
fig <- plot_ly(data = data.frame(reduced_P), x = ~v1, y = ~v2, type = "scatter")

fig


tag = rownames(reduced_P)
for(i in 1:(length(tag)-1)){
  tag[[i]] <- paste(strsplit(tag[[i]], "_")[[1]][-1],  collapse = " ")
}
tag

fig <- plot_ly(data = data.frame(reduced_P), x = ~v1, y = ~v2, type = "scatter", color = tag)

fig

# next thing to do :
# try to do pca while hiding some of the classes from the data :
# only plot breast & lung or any combinaison

# try to remove the most expressed genes from the heatmap to see if it become readable
#
#