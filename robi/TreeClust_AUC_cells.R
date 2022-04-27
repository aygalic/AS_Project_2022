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


row_names = reduced_data_treatment[,1]
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

M_ <- M

M_.e <- dist(M_, method='euclidean')
M_.m <- dist(M_, method='manhattan')
M_.c <- dist(M_, method='canberra')



# making every possible tree
M_.es <- hclust(M_.e, method='single')
M_.ea <- hclust(M_.e, method='average')
M_.ec <- hclust(M_.e, method='complete')

M_.ms <- hclust(M_.m, method='single')
M_.ma <- hclust(M_.m, method='average')
M_.mc <- hclust(M_.m, method='complete')

M_.cs <- hclust(M_.c, method='single')
M_.ca <- hclust(M_.c, method='average')
M_.cc <- hclust(M_.c, method='complete')





# We generate vectors of labels through the command cutree()
help(cutree)

# Fix k=2 clusters:
cluster.ec <- cutree(M_.ec, k=2) # euclidean-complete:
cluster.ec


library(factoextra) # clustering visualization

library(ggpubr) # ggarrange

# Making a silhouette plot with all metods
p1 <- fviz_nbclust(M_, FUN = hcut, method = "silhouette", hc_method ="complete", hc_metric="euclidian")
p2 <- fviz_nbclust(M_, FUN = hcut, method = "silhouette", hc_method ="complete", hc_metric="manhattan")
p3 <- fviz_nbclust(M_, FUN = hcut, method = "silhouette", hc_method ="complete", hc_metric="canberra")

p4 <- fviz_nbclust(M_, FUN = hcut, method = "silhouette", hc_method ="average", hc_metric="euclidian")
p5 <- fviz_nbclust(M_, FUN = hcut, method = "silhouette", hc_method ="average", hc_metric="manhattan")
p6 <- fviz_nbclust(M_, FUN = hcut, method = "silhouette", hc_method ="average", hc_metric="canberra")

p7 <- fviz_nbclust(M_, FUN = hcut, method = "silhouette", hc_method ="single", hc_metric="euclidian")
p8 <- fviz_nbclust(M_, FUN = hcut, method = "silhouette", hc_method ="single", hc_metric="manhattan")
p9 <- fviz_nbclust(M_, FUN = hcut, method = "silhouette", hc_method ="single", hc_metric="canberra")

silhouette_plot_all_method <- ggarrange(p1, p2, p3, p4, p5, p6, p7, p8 , p9,
                                        ncol = 3, nrow = 3)
x11()
silhouette_plot_all_method


# Making a ELBOW plot with all metods
e1 <- fviz_nbclust(M_, FUN = hcut, method = "wss", hc_method ="complete", hc_metric="euclidian")
e2 <- fviz_nbclust(M_, FUN = hcut, method = "wss", hc_method ="complete", hc_metric="manhattan")
e3 <- fviz_nbclust(M_, FUN = hcut, method = "wss", hc_method ="complete", hc_metric="canberra")

e4 <- fviz_nbclust(M_, FUN = hcut, method = "wss", hc_method ="average", hc_metric="euclidian")
e5 <- fviz_nbclust(M_, FUN = hcut, method = "wss", hc_method ="average", hc_metric="manhattan")
e6 <- fviz_nbclust(M_, FUN = hcut, method = "wss", hc_method ="average", hc_metric="canberra")

e7 <- fviz_nbclust(M_, FUN = hcut, method = "wss", hc_method ="single", hc_metric="euclidian")
e8 <- fviz_nbclust(M_, FUN = hcut, method = "wss", hc_method ="single", hc_metric="manhattan")
e9 <- fviz_nbclust(M_, FUN = hcut, method = "wss", hc_method ="single", hc_metric="canberra")

elbow_plot_all_method <- ggarrange(e1, e2, e3, e4, e5, e6, e7, e8 , e9,
                                   ncol = 3, nrow = 3)
x11()
elbow_plot_all_method

