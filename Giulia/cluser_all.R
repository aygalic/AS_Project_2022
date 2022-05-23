setwd("~/GitHub/AS_Project_2022")

#import data
data_patient= read.delim(file.path("Dataset/data_clinical_patient.txt"), header = TRUE, comment.char = '#')
data_treatment_auc= read.delim(file.path("Dataset", 'data_drug_treatment_auc.txt'), header = TRUE, comment.char = '#')
data_sample= read.delim(file.path("Dataset",'data_clinical_sample.txt'), header = TRUE, comment.char = '#')
original_data_mrna = read.delim(file.path("Dataset", "1_rpkm.txt"), header = TRUE, comment.char = '#')

cancer_data_treatment_auc<- data_treatment_auc
rownames(cancer_data_treatment_auc) <- data_treatment_auc$ENTITY_STABLE_ID
cancer_data_treatment_auc<- cancer_data_treatment_auc[,-1] #4 volte
#na

#col_names = reduced_data_treatment[1,]
perc_NA_col = colMeans(is.na(cancer_data_treatment_auc))*100

#row_names = reduced_data_treatment[,1]
perc_NA_row = rowMeans(is.na(cancer_data_treatment_auc))*100

mean(perc_NA_row) #20.85778

#I start with 50% tolerance on the rows
count=sum(perc_NA_row>50) #48

max_perc_row = max(perc_NA_row) #63.568
while(max_perc_row > 50)
{
  index = match(max_perc_row, perc_NA_row)
  cancer_data_treatment_auc = cancer_data_treatment_auc[-c(index),]
  perc_NA_row = perc_NA_row[-c(index)]
  max_perc_row = max(perc_NA_row)
}
#266 --> 218

#update the percentages
perc_NA_col = colMeans(is.na(cancer_data_treatment_auc))*100

#50% tolerance on the columns
max_perc_col = max(perc_NA_col)
while(max_perc_col > 50)
{
  index = match(max_perc_col, perc_NA_col)
  cancer_data_treatment_auc = cancer_data_treatment_auc[,-c(index)]
  perc_NA_col = perc_NA_col[-c(index)]
  max_perc_col = max(perc_NA_col)
}
#1065 --> 978

#update the percentages for rows
perc_NA_row = rowMeans(is.na(cancer_data_treatment_auc))*100

mean(perc_NA_row) #7.56271
mean(perc_NA_col) #7.56271
mean(is.na(cancer_data_treatment_auc))*100 #7.56271

#10% tolerance on the rows and I see if the mean percentage is 
max_perc_row = max(perc_NA_row)
while(max_perc_row > 10 && mean(perc_NA_row) > 2)
{
  index = match(max_perc_row, perc_NA_row)
  cancer_data_treatment_auc = cancer_data_treatment_auc[-c(index),]
  perc_NA_row = perc_NA_row[-c(index)]
  max_perc_row = max(perc_NA_row)
}

#218 --> 156
perc_NA_col = colMeans(is.na(cancer_data_treatment_auc))*100

mean(perc_NA_row) #4.108987
mean(perc_NA_col) #4.108987
mean(is.na(cancer_data_treatment_auc))*100 #4.108987


# max_perc_col = max(perc_NA_col)
# while(max_perc_col > 50)
# {
#   index = match(max_perc_col, perc_NA_col)
#   cancer_data_treatment_auc = cancer_data_treatment_auc[,-c(index)]
#   perc_NA_col = perc_NA_col[-c(index)]
#   max_perc_col = max(perc_NA_col)
# }
# #978-->974
# perc_NA_row = rowMeans(is.na(cancer_data_treatment_auc))*100
# 
# mean(perc_NA_row) #3.909
# mean(perc_NA_col) #3.909
# mean(is.na(cancer_data_treatment_auc))*100 #3.909
# 
# max_perc_col = max(perc_NA_col)
# while(max_perc_col > 40)
# {
#   index = match(max_perc_col, perc_NA_col)
#   cancer_data_treatment_auc = cancer_data_treatment_auc[,-c(index)]
#   perc_NA_col = perc_NA_col[-c(index)]
#   max_perc_col = max(perc_NA_col)
# }
# #978-->968
# perc_NA_row = rowMeans(is.na(cancer_data_treatment_auc))*100
# 
# mean(perc_NA_row) #3.672
# mean(perc_NA_col) #3.672
# mean(is.na(cancer_data_treatment_auc))*100 #3.672
sum(is.na(cancer_data_treatment_auc))
rows = rowMeans(cancer_data_treatment_auc, na.rm = TRUE)
sum(is.na(rows))
for(i in 1:156)
{
  for(j in 1:978)
  {
    if(is.na(cancer_data_treatment_auc[i,j])==TRUE)
    {
      cancer_data_treatment_auc[i,j] = rows[i]
    }
  }
}
sum(is.na(cancer_data_treatment_auc))



cancer_data_treatment_auc = t(cancer_data_treatment_auc)

M = as.data.frame(cancer_data_treatment_auc)
M <- as.matrix(cancer_data_treatment_auc)

de<- dist(M,method='euclidean')
dm<- dist(M,method='manhattan')
dc <- dist(M,method='canberra')

M_.es <- hclust(de, method='single')
M_.ea <- hclust(de, method='average')
M_.ec <- hclust(de, method='complete')

M_.ms <- hclust(dm, method='single')
M_.ma <- hclust(dm, method='average')
M_.mc <- hclust(dm, method='complete')

M_.cs <- hclust(dc, method='single')
M_.ca <- hclust(dc, method='average')
M_.cc <- hclust(dc, method='complete')


M_var= apply(M,2,var)
x11()
plot(M_var)
abline(h=0.06)

library(heatmaply)
M_0.06= M[,M_var>0.06]
heatmaply(data.matrix(M_0.06),hclustfun= hclust,dist_method='canberra' ,hclust_method = 'complete',k_row=7)




