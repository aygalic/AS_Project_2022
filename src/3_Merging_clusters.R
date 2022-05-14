setwd("~/OneDrive/polimi/COURSES/S8/APPLIED_STATS/AS_Project_2022")

source("src/utilities.R")

################################################
############### PREPARE THE DATA ###############
################################################

# steal Robi's work
M_AUC <- create_AUC_matrix()
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
reduced_M2_scaled <- create_reduced_mat(M2_)


# useful values
algos <- list(
  c("euclidian","single"),
  c("euclidian","average"),
  c("euclidian","complete"),
  c("manhattan","single"),
  c("manhattan","average"),
  c("manhattan","complete"),
  c("canberra","single"),
  c("canberra","average"),
  c("canberra","complete")
)
k = 10

# RPKM related stuff
v1.1 = reduced_M1_scaled$v1
v1.2 = reduced_M1_scaled$v2
names1 = rownames(reduced_M1_scaled)
size = length(v1.1)
k = 10

# AUC response related stuff
v2.1 = reduced_M2_scaled$v1
v2.2 = reduced_M2_scaled$v2
names2 = rownames(reduced_M2_scaled)
size = length(v2.1)


# FOR CONVEIGNANCE WE DECIDE TO PROJECT ON THE SAME REDUCED SPACE :
v2.1 = v1.1
v2.2 = v1.2


# FIRST WE BUILD THE FIRST TRACE 
# We une this one with kmeans since it's the only different algo

df1 <- data.frame(x = list(), y = list(), frame = list(), clust = list())
df2 <- data.frame(x = list(), y = list(), frame = list(), clust = list())


for(i in 1:k){
  df1 <- rbind(df1, data.frame(x = v1.1, y = v1.2, frame = rep(i, size), name = names1,
                               clust = kmeans(M1_, centers=i)$cluster))
  
  df2 <- rbind(df2, data.frame(x = v2.1, y = v2.2, frame = rep(i, size), name = names2,
                               clust = kmeans(M2_, centers=i)$cluster))
}





fig1 <- df1 %>% plot_ly()
fig2 <- df2 %>% plot_ly()

fig1 <- fig1 %>% add_markers(
  x = ~x, y = ~y, 
  hoverinfo = "text",
  text = ~name,
  frame = ~frame,
  color = ~clust,
  marker = list(colorscale = 'Viridis'),
  showlegend = F
) 

# since each figure has its own dataset there isn't a lot to change here
fig2 <- fig2 %>% add_markers(
  x = ~x, y = ~y, 
  hoverinfo = "text",
  text = ~name,
  frame = ~frame,
  color = ~clust,
  marker = list(colorscale = 'Viridis'),
  showlegend = F
) 

# THEN WE ADD A TRACE FOR EACH TREE ALGO
for(algo in algos){
  df1 <- data.frame( x = list(), y = list(), frame = list(), clust = list())
  df2 <- data.frame( x = list(), y = list(), frame = list(), clust = list())
  for(i in 1:k){
    df1 <- rbind(df1, data.frame(x = v1.1, y = v1.2, frame = rep(i ,size), name = names1,
                                 clust = hcut(M1_, i, hc_method = algo[2], hc_metric= algo[1])$cluster))
    
    df2 <- rbind(df2, data.frame(x = v2.1, y = v2.2, frame = rep(i ,size), name = names2,
                                 clust = hcut(M2_, i, hc_method = algo[2], hc_metric= algo[1])$cluster))
  }
  
  fig1 <- fig1 %>% add_markers(
    data = df1, x = ~x, y = ~y,
    text = ~name,
    hoverinfo = "text",
    frame = ~frame,
    color = ~clust,
    marker =list(colorscale = 'Viridis'),
    showlegend = F,
    visible = F
  )
  fig2 <- fig2 %>% add_markers(
    data = df2, x = ~x, y = ~y,
    text = ~name,
    hoverinfo = "text",
    frame = ~frame,
    color = ~clust,
    marker =list(colorscale = 'Viridis'),
    showlegend = F,
    visible = F
  )
}


fig <- subplot(fig1, fig2)
  
# add algorithm selection
fig <- fig %>% layout(
  title = "Comparing clustering Algo between the 2 datasets",
  xaxis = list(title = "PCA Axis 1"),
  yaxis = list(title = "PCA Axis 2"),
  updatemenus = list(
    list(
      y = 0.8,
      buttons = list(
        list(method = "restyle", args = list("visible", list(T, F, F, F, F, F, F, F, F, F)), label = "kmeans"),
        list(method = "restyle", args = list("visible", list(F, T, F, F, F, F, F, F, F, F)), label = "euclidian single"),       
        list(method = "restyle", args = list("visible", list(F, F, T, F, F, F, F, F, F, F)), label = "euclidian average"),        
        list(method = "restyle", args = list("visible", list(F, F, F, T, F, F, F, F, F, F)), label = "euclidian complete"),        
        list(method = "restyle", args = list("visible", list(F, F, F, F, T, F, F, F, F, F)), label = "manhattan single"),        
        list(method = "restyle", args = list("visible", list(F, F, F, F, F, T, F, F, F, F)), label = "manhattan average"),        
        list(method = "restyle", args = list("visible", list(F, F, F, F, F, F, T, F, F, F)), label = "manhattan complete"),        
        list(method = "restyle", args = list("visible", list(F, F, F, F, F, F, F, T, F, F)), label = "canberra single"),        
        list(method = "restyle", args = list("visible", list(F, F, F, F, F, F, F, F, T, F)), label = "canberra average"),        
        list(method = "restyle", args = list("visible", list(F, F, F, F, F, F, F, F, F, T)), label = "canberra complete")
      )
    )
  )
) 
# add the animation for the number of clusters
fig <- fig %>% animation_opts(0, easing = "elastic", redraw = TRUE) %>% hide_colorbar()
fig


saveWidget(fig, "output/aygalic/CLUSTERING_COMPARAISON.html", selfcontained = F, libdir = "lib")



