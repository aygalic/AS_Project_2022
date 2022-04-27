library(ggplot2)
library(plotly)  # interactive plots 
library(mvtnorm)
library(rgl)
library(car)
library(dplyr)


setwd("~/OneDrive/polimi/COURSES/S8/APPLIED_STATS/AS_Project_2022")


# steal Robi's work
source("robi/Kmeans_AUC_cells.R")
M1 <- M

# merge with mine
source("src/2_Kmeans_Clust_cells.R")
M2 <- M 

dev.off()


# select row from M1 that are present in M2
indexes1 <- rownames(M1) %in% rownames(M2)
M1_ <- M1 [indexes1,]
dim(M1_)

# select row from M2 that are present in M1
indexes2 <- rownames(M2) %in% rownames(M1)
M2_ <- M2[indexes2,]
dim(M2_)

# Do k_means on original matrix and keep only the labels of interest
result.AUC <- kmeans(M1, centers=k)$cluster[indexes1]
result.rpkm <- kmeans(M2, centers=k)$cluster[indexes2]

restult <- table(result.AUC, result.rpkm)

#make it correctly labbeled
rownames(restult) <- paste(rownames(restult),"_AUC")
colnames(restult) <- paste(rownames(restult),"_rpkm")
restult

#WHAT IS THE AMOUNT OF CORRECT LABELING ?
result.max <- apply(restult, 1, max)
result.sum <- apply(restult, 1, sum)
perf <- 100*result.max/result.sum
perf

cbind(restult, perf)

average_perf <- mean(perf)
average_perf











# now make a cool function out of that
contingency_table_avg_perf <- function(k = 2, fun = kmeans, 
                                     # param for htrees KEEP THEM NULL IF YOU USE KMEANS
                                     hc_method = NULL, hc_metric = NULL) 
  {
  
  result.AUC <- NULL
  result.rpkm <- NULL
  
  # case kmeans
  if(is.null(hc_method)){
    result.AUC <- fun(M1, centers=k)$cluster[indexes1]
    result.rpkm <- fun(M2, centers=k)$cluster[indexes2]
  }
  else{
    result.AUC <- fun(M1, k, hc_method = hc_method, hc_metric= hc_metric)$cluster[indexes1]
    result.rpkm <- fun(M2, k, hc_method = hc_method, hc_metric= hc_metric)$cluster[indexes2]
  }
  restult <- table(result.AUC, result.rpkm)
  
  #WHAT IS THE AMOUNT OF CORRECT LABELING ?
  result.max <- apply(restult, 1, max)
  result.sum <- apply(restult, 1, sum)

  return(mean(100*result.max/result.sum))
}

avg_perfs_kmeans <- c()
for(i in 1:10){
  avg_perfs_kmeans <- c(avg_perfs_kmeans, contingency_table_avg_perf(i, kmeans))
}
plot(avg_perfs_kmeans)


# this is a plot of how similar are the classifications between both datasets 
# depending on the number of cluster
avg_perfs_ce <- c()
avg_perfs_ae <- c()
avg_perfs_se <- c()

avg_perfs_cm <- c()
avg_perfs_am <- c()
avg_perfs_sm <- c()

avg_perfs_cc <- c()
avg_perfs_ac <- c()
avg_perfs_sc <- c()

j=10
for(i in 1:j){avg_perfs_ce <- c(avg_perfs_ce, contingency_table_avg_perf(i, hcut, hc_method = "complete", hc_metric ="euclidian"))}
for(i in 1:j){avg_perfs_ae <- c(avg_perfs_ae, contingency_table_avg_perf(i, hcut, hc_method = "average", hc_metric ="euclidian"))}
for(i in 1:j){avg_perfs_se <- c(avg_perfs_se, contingency_table_avg_perf(i, hcut, hc_method = "single", hc_metric ="euclidian"))}

for(i in 1:j){avg_perfs_cm <- c(avg_perfs_cm, contingency_table_avg_perf(i, hcut, hc_method = "complete", hc_metric ="manhattan"))}
for(i in 1:j){avg_perfs_am <- c(avg_perfs_am, contingency_table_avg_perf(i, hcut, hc_method = "average", hc_metric ="manhattan"))}
for(i in 1:j){avg_perfs_sm <- c(avg_perfs_sm, contingency_table_avg_perf(i, hcut, hc_method = "single", hc_metric ="manhattan"))}

for(i in 1:j){avg_perfs_cc <- c(avg_perfs_cc, contingency_table_avg_perf(i, hcut, hc_method = "complete", hc_metric ="canberra"))}
for(i in 1:j){avg_perfs_ac <- c(avg_perfs_ac, contingency_table_avg_perf(i, hcut, hc_method = "average", hc_metric ="canberra"))}
for(i in 1:j){avg_perfs_sc <- c(avg_perfs_sc, contingency_table_avg_perf(i, hcut, hc_method = "single", hc_metric ="canberra"))}


plot(avg_perfs_kmeans,type="l",col="black", ylim=c(0,100), 
     ylab = "AVG similar matching",
     xlab = "Number of clusters")

lines(avg_perfs_ce,col="indianred")
lines(avg_perfs_ae,col="indianred1")
lines(avg_perfs_se,col="indianred4")

lines(avg_perfs_cm,col="gold")
lines(avg_perfs_am,col="gold3")
lines(avg_perfs_sm,col="gold4")

lines(avg_perfs_cc,col="aquamarine")
lines(avg_perfs_ac,col="aquamarine3")
lines(avg_perfs_sc,col="aquamarine4")


legend("bottomleft", 95, 
       legend=c("kmeans",
                       "complete euclidian", "average euclidian", "single euclidian",
                       "complete manhattan", "average manhattan", "single manhattan",
                       "complete canberra", "average canberra", "single canberra"),
       col=c("black",
             "indianred", "indianred1","indianred4",
             "gold", "gold3", "gold4",
             "aquamarine", "aquamarine3", "aquamarine4"),
       lty=1:2, cex=0.8)
  

