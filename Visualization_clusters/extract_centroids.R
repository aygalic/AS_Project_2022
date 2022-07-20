load("breast_auc_data.Rdata")
load("selected_genes_data.Rdata")
breast_auc = t(breast_data_treatment_auc)
data_5 = na.aggregate(breast_auc) 

### MEANS
cell_means = apply(data_5,1,mean)

### PCA
cov2 = cov(data_5)
decomp2 = eigen(cov2)
P3 = as.matrix(decomp2$vectors[,1:3])
reduced_M_ = data_5%*%P3
colnames(reduced_M_) <- c("v1","v2", "v3")

data_plot = data.frame(reduced_M_)
data_plot$cell_means= cell_means
APER =1
while(APER>0.14){
  result.k <- kmeans(data_5, centers=3)
  cluster_group = result.k
  group <- as.factor(result.k$cluster)
  table(group)
  
  lda.fit <-  lda(cluster_group$cluster ~., data= data[,1:13])
  Lda.pred <-  predict(lda.fit , as.data.frame(data[,1:13]))
  
  n       <-   length(cluster_group$cluster)       # total number of observations
  ng      <-   table(cluster_group$cluster)       # number of obs. in each group
  group_names   <-   levels(group)      # name of groups
  g       <-   length(group_names)
  
  misc <- table(class.true=cluster_group$cluster, class.assigned=Lda.pred$class)
  
  #print(misc) #CONFUSION MATRIX
  
  APER=0
  for(gi in 1:g){
    APER <- APER + sum(misc[gi,-gi])/sum(misc[gi,]) * lda.fit$prior[gi]
  }
  
  print(APER)
  if(APER<0.10) break
}



# save(group, file="group_km.RData_3.RData")

# data_5
# names = c(NULL)
# for(i in 1:dim(result.k$centers)[1]){
#   center_data = result.k$centers[i,]
#   for(j in 1:dim(data_5)[1]){
#     if(norm(center_data - data_5[j,], type='2')<0.9){
#       names = c(names,rownames(data_5)[j])
#       break 
#       #print(rownames(data_5)[j])
#     }
#   }
# }
# names
# 
# mean_clusters
# 
# plot(2:dim(mean_clusters)[2],
#      mean_clusters[1,2:dim(mean_clusters)[2]],type="l", 
#      col="orange", pch = 19)
# 
# lines(2:dim(mean_clusters)[2],
#      mean_clusters[2,2:dim(mean_clusters)[2]],type="l", 
#      col="green")
# 
# lines(2:dim(mean_clusters)[2],
#       mean_clusters[3,2:dim(mean_clusters)[2]],type="l", 
#       col="red")
# legend("bottomleft", title="Mean of Group",
#        c("1","2","3"), fill=c("orange", "green", "red"), horiz=TRUE, cex=0.8)
# 
# mean_clusters[1,2:dim(mean_clusters)[2]]>mean_clusters[3,2:dim(mean_clusters)[2]]
# sum(mean_clusters[1,2:dim(mean_clusters)[2]]>mean_clusters[3,2:dim(mean_clusters)[2]])
# 

# -------------- DAY 19 JUL -----------
library(dplyr)
new_data <- to_return_2 %>% dplyr::select(-c('ZR7530_BREAST','BT549_BREAST'))
genes_data <- new_data[to_return_2$hugo_symbol %in% good_genes[1:13], ]
genes_hugo <- genes_data$hugo_symbol
data_pc <- t(genes_data %>% select(-c('hugo_symbol')))
colnames(data_pc) <- genes_hugo

# pca <- princomp(scale(data_pc), scores=T)
# P = pca$loading
#first component pca.std$loading[,1]
# k=3
# par(#mar = c(1,4,0,2),
#   mfrow = c(k,1))
# for(i in 1:k) 
#   barplot(P[,i], ylim = c(-1, 1))

k=3

colMeans(data_pc[group==2,])
colMeans(data_pc[group==3,])
par(#mar = c(1,4,0,2),
  mfrow = c(k,1))
for(i in 1:k) 
  barplot(colMeans(data_pc[group==i,]))

rbind(colMeans(data_pc[group==1,]),colMeans(data_pc[group==2,]), colMeans(data_pc[group==3,]))


library(tidyr)
df1 <- pivot_longer(as.data.frame(t(colMeans(data_pc[group==1,]))), cols=1:13, names_to = "Gene", values_to = "Expression")
df1$Group <- rep('1',13)

df2 <- pivot_longer(as.data.frame(t(colMeans(data_pc[group==2,]))), cols=1:13, names_to = "Gene", values_to = "Expression")
df2$Group <- rep('2',13)


df3 <- pivot_longer(as.data.frame(t(colMeans(data_pc[group==3,]))), cols=1:13, names_to = "Gene", values_to = "Expression")
df3$Group <- rep('3',13)

# df4 <- pivot_longer(as.data.frame(t(colMeans(data_pc[group==4,]))), cols=1:13, names_to = "Gene", values_to = "Expression")
# df4$Group <- rep('4',13)

total <- rbind(df1,df2,df3)

ggplot(total, aes(fill=Group, y=Expression, x=Gene)) + 
  geom_bar(position="dodge", stat="identity")
  