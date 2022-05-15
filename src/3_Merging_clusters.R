setwd("~/OneDrive/polimi/COURSES/S8/APPLIED_STATS/AS_Project_2022")

source("src/utilities.R")

################################################
############### PREPARE THE DATA ###############
################################################

# steal Robi's work
M_AUC <- create_AUC_matrix()
# TAKING ONLY FEW SAMPLES
M_AUC <- M_AUC[1:500,]
M1<- as.matrix(scale(M_AUC))

# merge with mine
M <- Build_matrix_for_multiple_cancer_types()$Mat
M_scaled <- as.matrix(scale(M))
M2 <- t(M_scaled)


# select row from M1 that are present in M2
indexes1 <- rownames(M1) %in% rownames(M2)
M1_ <- M1 [indexes1,]
dim(M1_)

# select row from M2 that are present in M1
indexes2 <- rownames(M2) %in% rownames(M1)
M2_ <- M2[indexes2,]
dim(M2_)





#################################################
################# ORIGINAL TEST #################
#################################################

# Do k_means on original matrix and keep only the labels of interest
k=10
result.AUC <- kmeans(M1, centers=k)$cluster[indexes1]
result.rpkm <- kmeans(M2, centers=k)$cluster[indexes2]

result <- table(result.AUC, result.rpkm)

#make it correctly labbeled
rownames(result) <- paste(rownames(result),"AUC")
colnames(result) <- paste(colnames(result),"rpkm")
result

#WHAT IS THE AMOUNT OF CORRECT LABELING ?
result.max <- apply(result, 1, max)
result.sum <- apply(result, 1, sum)
perf <- 100*result.max/result.sum
perf

cbind(result, perf)

average_perf <- mean(perf)
average_perf




# now make a cool function out of that
contingency_table <- function(k = 2, fun = kmeans, 
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
  result <- table(result.AUC, result.rpkm)
  
  #WHAT IS THE AMOUNT OF CORRECT LABELING ?
  result.max <- apply(result, 1, max)
  result.sum <- apply(result, 1, sum)


  perf <- 100*result.max/result.sum
  
  return(list(result = result, perf = perf, avg_perf=mean(perf)))
}

avg_perfs_kmeans <- c()
for(i in 1:10){
  avg_perfs_kmeans <- c(avg_perfs_kmeans, contingency_table(i, kmeans)$avg_perf)
}
plot(avg_perfs_kmeans)


# this is a plot of how similar are the classifications between both datasets 
# depending on the number of cluster
test_all_cluster_algo <- function(j=10){
  avg_perfs_kmeans <- c()
  
  avg_perfs_ce <- c()
  avg_perfs_ae <- c()
  avg_perfs_se <- c()
  
  avg_perfs_cm <- c()
  avg_perfs_am <- c()
  avg_perfs_sm <- c()
  
  avg_perfs_cc <- c()
  avg_perfs_ac <- c()
  avg_perfs_sc <- c()
  
  for(i in 1:j){avg_perfs_kmeans <- c(avg_perfs_kmeans, contingency_table(i, kmeans)$avg_perf)}
  
  for(i in 1:j){avg_perfs_ce <- c(avg_perfs_ce, contingency_table(i, hcut, hc_method = "complete", hc_metric ="euclidian")$avg_perf)}
  for(i in 1:j){avg_perfs_ae <- c(avg_perfs_ae, contingency_table(i, hcut, hc_method = "average", hc_metric ="euclidian")$avg_perf)}
  for(i in 1:j){avg_perfs_se <- c(avg_perfs_se, contingency_table(i, hcut, hc_method = "single", hc_metric ="euclidian")$avg_perf)}
  
  for(i in 1:j){avg_perfs_cm <- c(avg_perfs_cm, contingency_table(i, hcut, hc_method = "complete", hc_metric ="manhattan")$avg_perf)}
  for(i in 1:j){avg_perfs_am <- c(avg_perfs_am, contingency_table(i, hcut, hc_method = "average", hc_metric ="manhattan")$avg_perf)}
  for(i in 1:j){avg_perfs_sm <- c(avg_perfs_sm, contingency_table(i, hcut, hc_method = "single", hc_metric ="manhattan")$avg_perf)}
  
  for(i in 1:j){avg_perfs_cc <- c(avg_perfs_cc, contingency_table(i, hcut, hc_method = "complete", hc_metric ="canberra")$avg_perf)}
  for(i in 1:j){avg_perfs_ac <- c(avg_perfs_ac, contingency_table(i, hcut, hc_method = "average", hc_metric ="canberra")$avg_perf)}
  for(i in 1:j){avg_perfs_sc <- c(avg_perfs_sc, contingency_table(i, hcut, hc_method = "single", hc_metric ="canberra")$avg_perf)}
  
  
  
  
  n_clusters <- c(1:j)
  
  result <- data.frame(n_clusters, avg_perfs_kmeans,
                     avg_perfs_ce, avg_perfs_ae, avg_perfs_se,
                     avg_perfs_cm, avg_perfs_am, avg_perfs_sm,
                     avg_perfs_cc, avg_perfs_ac, avg_perfs_sc)
  
  
  fig <- plot_ly(result, x = ~n_clusters, y = ~avg_perfs_kmeans, name = 'kmeans', type = 'scatter', mode = 'lines+markers') 
  fig <- fig %>% add_trace(y = ~avg_perfs_ce, name = 'complete euclidian', mode = 'lines+markers') 
  fig <- fig %>% add_trace(y = ~avg_perfs_ae, name = 'average euclidian', mode = 'lines+markers') 
  fig <- fig %>% add_trace(y = ~avg_perfs_se, name = 'single euclidian', mode = 'lines+markers') 
  
  fig <- fig %>% add_trace(y = ~avg_perfs_cm, name = 'complete manhattan', mode = 'lines+markers') 
  fig <- fig %>% add_trace(y = ~avg_perfs_am, name = 'average manhattan', mode = 'lines+markers') 
  fig <- fig %>% add_trace(y = ~avg_perfs_sm, name = 'single manhattan', mode = 'lines+markers') 
  
  fig <- fig %>% add_trace(y = ~avg_perfs_cc, name = 'complete canberra"', mode = 'lines+markers') 
  fig <- fig %>% add_trace(y = ~avg_perfs_ac, name = 'average canberra', mode = 'lines+markers') 
  fig <- fig %>% add_trace(y = ~avg_perfs_sc, name = 'single canberra', mode = 'lines+markers') 
  
  return(list(result = result, fig = fig))
  

} 

result_10 <- test_all_cluster_algo(20)

result_10$fig

result_10$result








################################
######## Big brain plot ######## 
######## Big brain plot ######## 
######## Big brain plot ######## 
######## Big brain plot ######## 
######## Big brain plot ######## 
######## Big brain plot ######## 
################################

reduced_M1_scaled <- create_reduced_mat(M1_)
#reduced_M2_scaled <- create_reduced_mat(M2_)


# useful values
algos <- list(
  c("euclidian","single"), c("euclidian","average"), c("euclidian","complete"),
  c("manhattan","single"), c("manhattan","average"), c("manhattan","complete"),
  c("canberra", "single"), c("canberra", "average"), c("canberra", "complete")
)

all_algos <- list(
  "kmeans",
  "euclidian single",
  "euclidian average",
  "euclidian complete",
  "manhattan single",
  "manhattan average",
  "manhattan complete",
  "canberra single",
  "canberra average",
  "canberra complete"
)
n_algo = length(all_algos)
k = 10

# RPKM related stuff
v1 = reduced_M1_scaled$v1
v2 = reduced_M1_scaled$v2
names = rownames(reduced_M1_scaled)
size = length(v1)
k = 5

# We will use the projection axis of the RPKM dataset instead of the AUC one. 
# This shouldn't make a difference







# we prepare all the data
df1 <- data.frame(x = list(), y = list(), frame = list(), clust = list(), algo = list())
df2 <- data.frame(x = list(), y = list(), frame = list(), clust = list(), algo = list())


for(i in 1:k){
  df1 <- rbind(df1, data.frame(x = v1, y = v2, frame = rep(i, size), name = names,
                               clust = kmeans(M1_, centers=i)$cluster, 
                               algo = rep(all_algos[[1]], size)))
  
  df2 <- rbind(df2, data.frame(x = v1, y = v2, frame = rep(i, size), name = names,
                               clust = kmeans(M2_, centers=i)$cluster, 
                               algo = rep(all_algos[[1]], size)))
}


for(j in 2:n_algo){
  for(i in 1:k){
    df1 <- rbind(df1, data.frame(x = v1, y = v2, frame = rep(i ,size), name = names,
                 clust = hcut(M1_, i, hc_method = algos[[j-1]][2], hc_metric= algos[[j-1]][1])$cluster, 
                 algo = rep(all_algos[[j]], size)))
    
    df2 <- rbind(df2, data.frame(x = v1, y = v2, frame = rep(i ,size), name = names,
                 clust = hcut(M2_, i, hc_method = algos[[j-1]][2], hc_metric= algos[[j-1]][1])$cluster,
                 algo = rep(all_algos[[j]], size)))
  }
}




# we build the plots
fig1 <- plot_ly()
fig2 <- plot_ly()
for(alg in all_algos){
  # optimization
  VISIBLE = ifelse(alg=="kmeans", T, F)
  df1_ <- df1[df1$algo == alg,]
  df2_ <- df2[df2$algo == alg,]
  #dubug
  print(dim(df1_))
  print(alg)
  print(VISIBLE)
  
  fig1 <- fig1 %>% add_markers(
    data = df1_,
    x = ~x, y = ~y,
    text = paste(df1_$name, alg),
    hoverinfo = "text",
    frame = ~frame,
    color = ~clust,
    marker =list(colorscale = 'Picnic'),
    #showlegend = F, 
    visible = VISIBLE
  )
  fig2 <- fig2 %>% add_markers(
    data = df2_,
    x = ~x, y = ~y,
    text = paste(df2_$name, alg),
    hoverinfo = "text",
    frame = ~frame,
    color = ~clust,
    marker =list(colorscale = 'Picnic'),
    #showlegend = F,
    visible = VISIBLE
  )
}
fig1






METHOD = "update"
BTN = list(
  list(method = METHOD, args = list(list(visible = c(T, F, F, F, F, F, F, F, F, F)), list(title = "1")), label = "kmeans"),
  list(method = METHOD, args = list(list(visible = c(F, T, F, F, F, F, F, F, F, F)), list(title = "2")), label = "euclidian single"),       
  list(method = METHOD, args = list(list(visible = c(F, F, T, F, F, F, F, F, F, F)), list(title = "3")), label = "euclidian average"),        
  list(method = METHOD, args = list(list(visible = c(F, F, F, T, F, F, F, F, F, F)), list(title = "4")), label = "euclidian complete"),        
  list(method = METHOD, args = list(list(visible = c(F, F, F, F, T, F, F, F, F, F)), list(title = "5")), label = "manhattan single"),        
  list(method = METHOD, args = list(list(visible = c(F, F, F, F, F, T, F, F, F, F)), list(title = "6")), label = "manhattan average"),        
  list(method = METHOD, args = list(list(visible = c(F, F, F, F, F, F, T, F, F, F)), list(title = "7")), label = "manhattan complete"),        
  list(method = METHOD, args = list(list(visible = c(F, F, F, F, F, F, F, T, F, F)), list(title = "8")), label = "canberra single"),        
  list(method = METHOD, args = list(list(visible = c(F, F, F, F, F, F, F, F, T, F)), list(title = "9")), label = "canberra average"),        
  list(method = METHOD, args = list(list(visible = c(F, F, F, F, F, F, F, F, F, T)), list(title = "10")), label = "canberra complete")
)

BTN1 = list(
  list(method = METHOD, args = list(list(visible = c(T, F, F, F, F, F, F, F, F, F, F, T, F, F, F, F, F, F, F, F, F, F)), list(title = "1")), label = "kmeans"),
  list(method = METHOD, args = list(list(visible = c(F, T, F, F, F, F, F, F, F, F, F, F, T, F, F, F, F, F, F, F, F, F)), list(title = "2")), label = "euclidian single"),       
  list(method = METHOD, args = list(list(visible = c(F, F, T, F, F, F, F, F, F, F, F, F, F, T, F, F, F, F, F, F, F, F)), list(title = "3")), label = "euclidian average"),        
  list(method = METHOD, args = list(list(visible = c(F, F, F, T, F, F, F, F, F, F, F, F, F, F, T, F, F, F, F, F, F, F)), list(title = "4")), label = "euclidian complete"),        
  list(method = METHOD, args = list(list(visible = c(F, F, F, F, T, F, F, F, F, F, F, F, F, F, F, T, F, F, F, F, F, F)), list(title = "5")), label = "manhattan single"),        
  list(method = METHOD, args = list(list(visible = c(F, F, F, F, F, T, F, F, F, F, F, F, F, F, F, F, T, F, F, F, F, F)), list(title = "6")), label = "manhattan average"),        
  list(method = METHOD, args = list(list(visible = c(F, F, F, F, F, F, T, F, F, F, F, F, F, F, F, F, F, T, F, F, F, F)), list(title = "7")), label = "manhattan complete"),        
  list(method = METHOD, args = list(list(visible = c(F, F, F, F, F, F, F, T, F, F, F, F, F, F, F, F, F, F, T, F, F, F)), list(title = "8")), label = "canberra single"),        
  list(method = METHOD, args = list(list(visible = c(F, F, F, F, F, F, F, F, T, F, F, F, F, F, F, F, F, F, F, T, F, F)), list(title = "9")), label = "canberra average"),        
  list(method = METHOD, args = list(list(visible = c(F, F, F, F, F, F, F, F, F, T, F, F, F, F, F, F, F, F, F, F, T, F)), list(title = "10")), label = "canberra complete")
)


# add algorithm selection

fig <- subplot(fig1, fig2) %>% animation_opts(0, easing = "linear", redraw = F) 
fig <- fig %>% layout(
  title = "Comparing clustering Algo between the 2 datasets",
  xaxis = list(title = "PCA Axis 1"),
  yaxis = list(title = "PCA Axis 2"),
  updatemenus = list(
    list(
      y = 0.8,
      active = 0,
      type= 'buttons',
      buttons = BTN1
    )
  )
) 


#fig %>% hide_colorbar()
fig

saveWidget(fig, "output/aygalic/CLUSTERING_COMPARAISON.html", selfcontained = F, libdir = "lib")



