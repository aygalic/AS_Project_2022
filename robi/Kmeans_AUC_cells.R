library(ggplot2)
library(plotly)  
library(mvtnorm)
library(rgl)
library(car)
library(factoextra)
library(NbClust)


#import data
data_patient = read.delim(file.path("Dataset", "data_clinical_patient.txt"), header = TRUE, comment.char = '#')

data_treatment_auc = read.delim(file.path("Dataset", 'data_drug_treatment_auc.txt'), header = TRUE, comment.char = '#')

data_sample = read.delim(file.path("Dataset",'data_clinical_sample.txt'), header = TRUE, comment.char = '#')

sample_ID = data_sample$SAMPLE_ID

indexes = match(sample_ID, colnames(data_treatment_auc))
indexes

sum(is.na(indexes)) 

cancer_data_treatment_auc = data_treatment_auc[, c(na.omit(indexes)) ]
rownames(cancer_data_treatment_auc) = data_treatment_auc$ENTITY_STABLE_ID

dim(cancer_data_treatment_auc)
#we have 266 treatments and 1063 cell lines

#remove NA
mean(rowMeans(is.na(cancer_data_treatment_auc))*100) #20.85% of NA

#Remove NA
data_treatment_auc = cancer_data_treatment_auc
col_names = data_treatment_auc[1,]
perc_NA_col = colMeans(is.na(data_treatment_auc))*100


row_names = data_treatment_auc[,1]
perc_NA_row = rowMeans(is.na(data_treatment_auc))*100

sum(perc_NA_row>50)
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

mean(perc_NA_row)
mean(perc_NA_col)

sum(perc_NA_row>10)

max_perc_row = max(perc_NA_row)
while(max_perc_row > 8)
{
  index = match(max_perc_row, perc_NA_row)
  data_treatment_auc = data_treatment_auc[-c(index),]
  perc_NA_row = perc_NA_row[-c(index)]
  max_perc_row = max(perc_NA_row)
}

mean(perc_NA_row)
mean(perc_NA_col)
mean(is.na(data_treatment_auc))*100 # 3.466

dim(data_treatment_auc) #137 rows, 976 columns

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
View(data_treatment_auc)

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

View(M)

M <- as.matrix(scale(M))
result.k <- kmeans(M, centers=3)

names(result.k)

result.k$cluster      # labels of clusters
result.k$centers      # centers of the clusters
result.k$totss        # tot. sum of squares
result.k$withinss     # sum of squares within clusters
result.k$tot.withinss # sum(sum of squares nei cluster)
result.k$betweenss    # sum of squares between clusters
result.k$size         # dimention of the clusters

#plot on the first 2 treatments
x11()
plot(M, col = result.k$cluster+1)

#plot on the first 3 treatments
open3d()
plot3d(M, size=3, col=result.k$cluster+1, aspect = F) 
points3d(result.k$centers,size=10)

# k = 4
result.k4 <- kmeans(M, centers=4)
x11()
plot(M, col = result.k4$cluster+1)

#find the best number of clusters
fviz_nbclust(M, FUN = kmeans, method = "silhouette") 
fviz_nbclust(M, FUN = kmeans, method = "wss")
#3 minimum optimal number of clusters, up to 7 clusters maximum


