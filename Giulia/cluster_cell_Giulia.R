#reduced_data_treatment is obtained running the code analysis_Giulia until line 132
#It's the dataset with only breast cancer and without NA
cancer_data_treatment_auc = t(reduced_data_treatment)

M = as.data.frame(cancer_data_treatment_auc)
M <- as.matrix(cancer_data_treatment_auc)

#dissimilarity
de<- dist(M,method='euclidean')
dm<- dist(M,method='manhattan')
dc <- dist(M,method='canberra')

x11()
par(mfrow=c(1,3))
image(as.matrix(de),main=' euclidean')
image(as.matrix(dm), main='manhattan')
image(as.matrix(dc),main='canberra')

#hclust 
M_.es <- hclust(de, method='single')
M_.ea <- hclust(de, method='average')
M_.ec <- hclust(de, method='complete')

M_.ms <- hclust(dm, method='single')
M_.ma <- hclust(dm, method='average')
M_.mc <- hclust(dm, method='complete')

M_.cs <- hclust(dc, method='single')
M_.ca <- hclust(dc, method='average')
M_.cc <- hclust(dc, method='complete')

#dendogram
x11()
par(mfrow=c(3,3))
plot(M_.es, main='euclidean-single', hang=-0.1, xlab='', labels=F, cex=0.6, sub='')
plot(M_.ec, main='euclidean-complete', hang=-0.1, xlab='', labels=F, cex=0.6, sub='')
plot(M_.ea, main='euclidean-average', hang=-0.1, xlab='', labels=F, cex=0.6, sub='')

plot(M_.ms, main='euclidean-single', hang=-0.1, xlab='', labels=F, cex=0.6, sub='')
plot(M_.mc, main='euclidean-complete', hang=-0.1, xlab='', labels=F, cex=0.6, sub='')
plot(M_.ma, main='euclidean-average', hang=-0.1, xlab='', labels=F, cex=0.6, sub='')

plot(M_.cs, main='euclidean-single', hang=-0.1, xlab='', labels=F, cex=0.6, sub='')
plot(M_.cc, main='euclidean-complete', hang=-0.1, xlab='', labels=F, cex=0.6, sub='')
plot(M_.ca, main='euclidean-average', hang=-0.1, xlab='', labels=F, cex=0.6, sub='')



#set the number of clusters
k=3

cluster.ec <- cutree(M_.ec, k=k) # euclidean-complete
cluster.es <- cutree(M_.es, k=k) # euclidean-single
cluster.ea <- cutree(M_.ea, k=k) # euclidean-average

cluster.ms <- cutree(M_.ms, k=k) # manhattan-single
cluster.ma <- cutree(M_.ma, k=k) # manhattan-average
cluster.mc <- cutree(M_.mc, k=k) # manhattan-complete

cluster.cs <- cutree(M_.cs, k=k) # canberra-single
cluster.ca <- cutree(M_.ca, k=k) # canberra-average
cluster.cc <- cutree(M_.cc, k=k) # canberra-complete

#table ...
table(label.true = rownames(M), label.cluster = cluster.es)

x11()
par(mfrow=c(3,3))
plot(M, main="euclidean-single", col=cluster.es, pch=19)
plot(M, main="euclidean-average", col=cluster.ea, pch=19)
plot(M, main="euclidean-complete", col=cluster.ec, pch=19)

plot(M, main="manhattan-single", col=cluster.ms, pch=19)
plot(M, main="manhattan-average", col=cluster.ma, pch=19)
plot(M, main="manhattan-complete", col=cluster.mc, pch=19)

plot(M, main="canberra-single", col=cluster.cs, pch=19)
plot(M, main="canberra-average", col=cluster.ca, pch=19)
plot(M, main="canberra-complete",col=cluster.cc, pch=19)

#plot on principal axes
#Aygalic function
reduced_M_scaled = create_reduced_mat(M,2)

x11()
par(mfrow=c(3,3))
plot(data.frame(reduced_M_scaled), main = 'euclidean single',   col=cluster.es, pch=16, asp=1)
plot(data.frame(reduced_M_scaled), main = 'euclidean average',  col=cluster.ea, pch=16, asp=1)
plot(data.frame(reduced_M_scaled), main = 'euclidean complete', col=cluster.ec, pch=16, asp=1)
plot(data.frame(reduced_M_scaled), main = 'manhattan single',   col=cluster.ms, pch=16, asp=1)
plot(data.frame(reduced_M_scaled), main = 'manhattan average',  col=cluster.ma, pch=16, asp=1)
plot(data.frame(reduced_M_scaled), main = 'manhattan complete', col=cluster.mc, pch=16, asp=1)
plot(data.frame(reduced_M_scaled), main = 'canberra single',    col=cluster.cs, pch=16, asp=1)
plot(data.frame(reduced_M_scaled), main = 'canberra average',   col=cluster.ca, pch=16, asp=1)
plot(data.frame(reduced_M_scaled), main = 'canberra complete',  col=cluster.cc, pch=16, asp=1)

dev.off()


#Kmeans

result.k <- kmeans(M, centers=3)
 
x11()
par(mfrow=c(1,2))
plot(M, col = result.k$cluster+1)
plot(reduced_M_scaled, col=result.k$cluster+1)



#compare
coph.es <- cophenetic(M_.es)
coph.ec <- cophenetic(M_.ec)
coph.ea <- cophenetic(M_.ea)

coph.ms <- cophenetic(M_.ms)
coph.mc <- cophenetic(M_.mc)
coph.ma <- cophenetic(M_.ma)

coph.cs <- cophenetic(M_.cs)
coph.cc <- cophenetic(M_.cc)
coph.ca <- cophenetic(M_.ca)

x11()
layout(rbind(c(0,1,0),c(2,3,4)))
image(as.matrix(de))
image(as.matrix(coph.es))
image(as.matrix(coph.ec))
image(as.matrix(coph.ea))

x11()
layout(rbind(c(0,1,0),c(2,3,4)))
image(as.matrix(dm))
image(as.matrix(coph.ms))
image(as.matrix(coph.mc))
image(as.matrix(coph.ma))

x11()
layout(rbind(c(0,1,0),c(2,3,4)))
image(as.matrix(dc))
image(as.matrix(coph.cs))
image(as.matrix(coph.cc))
image(as.matrix(coph.ca))


#I use kmeans from now

#Select the cell line of the clusters

groups <- result.k$cluster
matrix_group= data.frame(groups)

#index in the AUC matrix
index_1= which(matrix_group==1)
index_2= which(matrix_group==2)
index_3= which(matrix_group==3)
#index_4= which(matrix_group==4)
#index_5= which(matrix_group==5)

name=rownames(matrix_group)
name1=name[index_1]
name2=name[index_2]
name3=name[index_3]
#name4=name[index_4]
#name5=name[index_5]

M1=M[index_1,]
M2=M[index_2,]
M3=M[index_3,]
M_new= rbind(M1,M2,M3)
library(heatmaply)
heatmaply(data.matrix(M_new))

#variance
M_var= apply(M,2,var)
x11()
plot(M_var)

M_0.06= M[,M_var>0.06]
heatmaply(data.matrix(M_0.06))

#rpkm dataset
R=as.data.frame(mrna_data)
R= t(R)
R=scale(R)

sum(is.na(R))
for(i in 1:length(colnames(R)))
  for(j in 1:length(rownames(R)))
    if(is.na(R[j,i]))
      R=R[,-c(i)]

rpkm_col= rownames(R)
rpkm_1<-NULL
rpkm_2 <-NULL
rpkm_3 <-NULL
for ( i in 1:22)
  rpkm_1= c(rpkm_1,which(rpkm_col==name1[i]))
# rpkm_1 are 12 --> one missing
for( i in 1:18)
{
  rpkm_2= c(rpkm_2,which(rpkm_col==name2[i]))
}
# 22
for ( i in 1:3)
{
  rpkm_3= c(rpkm_3,which(rpkm_col==name3[i]))
}
# 8


reduced_R_scaled= create_reduced_mat(R,2)

#plot the cluster get from AUC on the reduced_R_scaled axes
X11()
par(mfrow=c(1,2))
plot(R, main= 'rpkm',col=groups)
plot(reduced_R_scaled, main='rpkm on PC', col=groups)

#expression
#scaling
#most expressed genes

index_match_1= match(name1, rownames(R))
index_match_2= match(name2, rownames(R))
index_match_3= match(name3, rownames(R))

#rpkm matrix of the 3 clusters
R1 = R[index_match_1,]
R2=  R[index_match_2,]
R3=  R[index_match_3,] 

R_new= rbind(R1,R2,R3)
library(heatmaply)
heatmaply(data.matrix(R_new))

#most expressed genes in R1
col1<-NULL
for(i in 1:length(index_match_1))
{
  for(j in 1:length(colnames(R1)))
    if(abs(R1[i,j])>0.7)
      col1[j]=1
  else
    col1[j]=2
}
reduced_R1= create_reduced_mat(R1,2)
X11()
par(mfrow=c(1,2))
plot(R1,col=col1)
plot(reduced_M_scaled, col=col1)
max(R1) #7.21225
min(R1) #-2.2371
sum(R1>0) #6115
sum(R1<0) #13091
sum(R1>5) #72
sum(R1<5) #19134
R1_mean=apply(R1,2,mean)
R1_var= apply(R1,2,var)
x11()
par(mfrow=c(1,2))
plot(R1_var)
plot(R1_mean)

heatmaply(data.matrix(R1))


col2<-NULL
for(i in 1:length(index_match_2))
{
  for(j in 1:length(colnames(R2)))
    if(abs(R2[i,j])>0.5)
      col2[j]=1
  else
    col2[j]=2
}
reduced_R2= create_reduced_mat(R2,2)
X11()
par(mfrow=c(1,2))
plot(R2,col=col2)
plot(reduced_R2,col=col2)
max(R2) #7.212
min(R2) #-2.2846

col3<-NULL
for(i in 1:length(index_match_3))
{
  for(j in 1:length(colnames(R3)))
    if(abs(R3[i,j])>0.8)
      col2[j]=1
  else
    col3[j]=2
}
reduced_R3= create_reduced_mat(R3,2)
X11()
par(mfrow=c(1,2))
plot(R3,col=col3)
plot(reduced_R3,col=col3)

max(R3) #7.0511
min(R3) #-1.723492




library(fpc)

db <- dbscan(M,eps=5,MinPts = 2)
db$cluster
db

