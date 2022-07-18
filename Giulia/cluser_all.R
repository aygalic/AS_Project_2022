setwd("~/GitHub/AS_Project_2022")

#import data (line 3-116, only RUN)
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

#Plot of the variance of the treatment, i take as  treshold 0.06, feel free to change
M_var= apply(M,2,var)
x11()
plot(M_var)
abline(h=0.06)
M_0.06= M[,M_var>0.06]

#All cluster 
de_0.06<- dist(M_0.06,method='euclidean')
dm_0.06<- dist(M_0.06,method='manhattan')
dc_0.06 <- dist(M_0.06,method='canberra')

M_.es_0.06 <- hclust(de, method='single')
M_.ea_0.06 <- hclust(de, method='average')
M_.ec_0.06 <- hclust(de, method='complete')

M_.ms_0.06 <- hclust(dm, method='single')
M_.ma_0.06 <- hclust(dm, method='average')
M_.mc_0.06 <- hclust(dm, method='complete')

M_.cs_0.06 <- hclust(dc, method='single')
M_.ca_0.06 <- hclust(dc, method='average')
M_.cc_0.06 <- hclust(dc, method='complete')

#dendogram
x11()
par(mfrow=c(1,3))
plot(M_.es_0.06,labels= FALSE,main="ES")
rect.hclust(M_.es_0.06, k=7)
plot(M_.ea_0.06,labels= FALSE,main="EA")
rect.hclust(M_.ea_0.06, k=7)
plot(M_.ec_0.06,labels= FALSE,main="EC")
rect.hclust(M_.ec_0.06, k=7)


x11()
par(mfrow=c(1,3))
plot(M_.ms_0.06,labels= FALSE,main="MS")
rect.hclust(M_.ms_0.06, k=7)
plot(M_.ma_0.06,labels= FALSE,main="MA")
rect.hclust(M_.ma_0.06, k=7)
plot(M_.mc_0.06,labels= FALSE,main="MC")
rect.hclust(M_.mc_0.06, k=7)

x11()
par(mfrow=c(1,3))
plot(M_.cs_0.06,labels= FALSE,main="CS")
rect.hclust(M_.cs_0.06, k=7)
plot(M_.ca_0.06,labels= FALSE,main="CA")
rect.hclust(M_.ca_0.06, k=7)
plot(M_.cc_0.06,labels= FALSE,main="CC")
rect.hclust(M_.cc_0.06, k=7)

#single doesn't make sense to me

#compare the average
x11()
par(mfrow=c(1,3))
plot(M_.ea_0.06,labels= FALSE,main="EA")
rect.hclust(M_.ea_0.06, k=7)
plot(M_.ma_0.06,labels= FALSE)
rect.hclust(M_.ma_0.06, k=7)
plot(M_.ca_0.06,labels= FALSE)
rect.hclust(M_.ca_0.06, k=7)

#compare the complete
x11()
par(mfrow=c(1,3))
plot(M_.ec_0.06,labels= FALSE,main="EC")
rect.hclust(M_.ec_0.06, k=7)
plot(M_.mc_0.06,labels= FALSE,main="MC")
rect.hclust(M_.mc_0.06, k=7)
plot(M_.cc_0.06,labels= FALSE,main="CC")
rect.hclust(M_.cc_0.06, k=7)


k=7

cluster.ec <- cutree(M_.ec_0.06, k=k) # euclidean-complete
cluster.es <- cutree(M_.es_0.06, k=k) # euclidean-single
cluster.ea <- cutree(M_.ea_0.06, k=k) # euclidean-average

cluster.ms <- cutree(M_.ms_0.06, k=k) # manhattan-single
cluster.ma <- cutree(M_.ma_0.06, k=k) # manhattan-average
cluster.mc <- cutree(M_.mc_0.06, k=k) # manhattan-complete

cluster.cs <- cutree(M_.cs_0.06, k=k) # canberra-single
cluster.ca <- cutree(M_.ca_0.06, k=k) # canberra-average
cluster.cc <- cutree(M_.cc_0.06, k=k) # canberra-complete

#I take the canberra-complete
x11()
plot(cluster.cc)

library(heatmaply)
heatmaply(data.matrix(M_0.06),hclustfun= hclust,dist_method='canberra' ,hclust_method = 'complete',k_row=7)

groups <- cluster.cc
matrix_group= data.frame(cluster.cc)

#index in the AUC matrix
index_1= which(matrix_group==1)
index_2= which(matrix_group==2)
index_3= which(matrix_group==3)
index_4= which(matrix_group==4)
index_5= which(matrix_group==5)
index_6= which(matrix_group==6)
index_7= which(matrix_group==7)

name=rownames(matrix_group)
name1=name[index_1]
name2=name[index_2]
name3=name[index_3]
name4=name[index_4]
name5=name[index_5]
name6=name[index_6]
name7=name[index_7]

#Understand if in the group we have the same part of the body or if they are different
#In group 7 we have only HAEMATOPOIETIC_AND_LYMPHOID_TISSUE
#In group 6 we have HAEMATOPOIETIC_AND_LYMPHOID_TISSUE, one OVARY, two LUNG, two BONE

M7=M[index_7,]
heatmap(M7)
