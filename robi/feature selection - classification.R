breast_rpkm_data = read.delim(file.path("Dataset", "2_rpkm.txt"), header = TRUE, comment.char = '#')
breast_rpkm_data = t(breast_rpkm_data)
#colnames(breast_rpkm_data) = breast_rpkm_data[1,]
#breast_rpkm_data = breast_rpkm_data[-1,]

source("src/utilities.R")
source("robi/analysis_robi.rmd")

breast_auc_data = cancer_data_treatment_auc

#set the mean as value fo NAs
rows = rowMeans(breast_auc_data, na.rm = TRUE)

for(i in 1:142)
{
  for(j in 1:43)
  {
    if(is.na(breast_auc_data[i,j])==TRUE)
    {
      breast_auc_data[i,j] = rows[i]
    }
  }
}
sum(is.na(breast_auc_data))
#no more NAs in the dataframe

#transpose to get cell lines in the rows
breast_auc_data = t(breast_auc_data)

M1 = as.matrix(breast_auc_data)
M2 = as.matrix(breast_rpkm_data)

# select row from M1 that are present in M2
indexes1 <- rownames(M1) %in% rownames(M2)
M1_ <- M1[indexes1,]
dim(M1_)

# select row from M2 that are present in M1
indexes2 <- rownames(M2) %in% rownames(M1)
M2_ <- M2[indexes2,]
dim(M2_)

result.k <- kmeans(M1_, centers=3)

reduced_M1_scaled = create_reduced_mat(M1_,2)

x11()
par(mfrow=c(1,2))
plot(M1_, col = result.k$cluster+1) #only the first 2 treatments
plot(reduced_M1_scaled, col=result.k$cluster+1, pch = 16, xlab = "PC1", ylab = "PC2") #on PC1/PC2 space

true_labels <- result.k$cluster

label_1 = which(true_labels == 1)
label_2 = which(true_labels == 2)
label_3 = which(true_labels == 3)


