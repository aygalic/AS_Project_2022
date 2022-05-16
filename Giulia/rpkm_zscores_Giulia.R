#rpkm Z-SCORES

data_patient_= read.delim(file.path("Dataset", "data_clinical_patient.txt"), header = TRUE, comment.char = '#')
data_sample= read.delim(file.path("Dataset",'data_clinical_sample.txt'), header = TRUE, comment.char = '#')
original_data_mrna = read.delim(file.path("Dataset", "data_mrna_seq_rpkm.txt"), header = TRUE, comment.char = '#', nrows=100)

rownames(original_data_mrna) <- original_data_mrna$Hugo_Symbol
original_data_mrna[,1] <-NULL

#select cells from breast carcinoma patients
cancer_name= "Invasive Breast Carcinoma"
selected_values = data_sample[data_sample$"CANCER_TYPE_DETAILED"==cancer_name,]
selected_cells = selected_values$SAMPLE_ID
selected_cells = na.omit(selected_cells)
indexes_mrna = match(selected_cells, colnames(original_data_mrna))
mrna_data = original_data_mrna[, c(1, na.omit(indexes_mrna)) ]

#duplicate

sum(duplicated(rownames(mrna_data)))

#NA value: gene without diploid sample
sum(is.na(mrna_data)) #162 --> 3 rows completely NA
null_row <- 1:3
j=1
for (i in 1:100)
  if (is.na(mrna_data[i,1])==TRUE)
  {
    null_row[j]=i
    j= j+1
  }
mrna_data= mrna_data[-c(null_row),]


#heatmap
library(heatmaply)
heatmaply(mrna_data)

#how much above 10?
sum(mrna_data>10) #3; I put this three equal to 10
for (i in 1:dim(mrna_data)[1])
  for ( j in 1:dim(mrna_data)[2])
    if (mrna_data[i,j]>10)
      mrna_data[i,j]=10
#I'm not sure about this part but we have only 3 elements of the matrix really high
#We have a scale between 0 and 20 and only 3 are above ten

#PCA

pc <- princomp(mrna_data,scores=T)
summary(pc)

#sdv
dev.new()
plot(pc, las=2, main='Principal components', ylim=c(0,1))
barplot(sapply(mrna_data,sd)^2, las=2, main='Original Variables', ylim=c(0,0.1), ylab='Variances',col='pink')
plot(cumsum(pc$sd^2)/sum(pc$sd^2), type='b', xlab='number of components', 
     ylab='contribution to the total variance', ylim=c(0,1))
abline(h=0.45, col='red') #12 components for 80%
                          #3 components for 45%
dev.off()


#loadings and treatment
load.treatment <- pc$loadings
scores.treatment <- pc$scores

load.treatment[load.treatment>0.5]
dev.new()
for(i in 1:3)barplot(load.treatment[,i], ylim = c(-0.5, 0.5), main=paste('Loadings PC ',i,sep=''))
for(i in 1:3)barplot(scores.treatment[,i], ylim = c(-15, 10), main=paste('Scores PC ',i,sep=''))

dev.off()

library(pca3d)
pca3d(pc)

pca2d(pc)

#Find a way to color the pc if it is possible

