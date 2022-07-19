data_5
names = c(NULL)
for(i in 1:dim(result.k$centers)[1]){
  center_data = result.k$centers[i,]
  for(j in 1:dim(data_5)[1]){
    if(norm(center_data - data_5[j,], type='2')<0.9){
      names = c(names,rownames(data_5)[j])
      break 
      #print(rownames(data_5)[j])
    }
  }
}
names


plot(2:dim(mean_clusters)[2],
     mean_clusters[1,2:dim(mean_clusters)[2]],type="l", 
     col="orange", pch = 19)

lines(2:dim(mean_clusters)[2],
     mean_clusters[2,2:dim(mean_clusters)[2]],type="l", 
     col="green")

lines(2:dim(mean_clusters)[2],
      mean_clusters[3,2:dim(mean_clusters)[2]],type="l", 
      col="red")
legend("bottomleft", title="Mean of Group",
       c("1","2","3"), fill=c("orange", "green", "red"), horiz=TRUE, cex=0.8)

mean_clusters[1,2:dim(mean_clusters)[2]]>mean_clusters[3,2:dim(mean_clusters)[2]]
sum(mean_clusters[1,2:dim(mean_clusters)[2]]>mean_clusters[3,2:dim(mean_clusters)[2]])
  