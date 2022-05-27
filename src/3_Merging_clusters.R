setwd("~/OneDrive/polimi/COURSES/S8/APPLIED_STATS/AS_Project_2022")

source("src/utilities.R")

################################################
############### PREPARE THE DATA ###############
################################################

# steal Robi's work
M_AUC <- create_AUC_matrix()
# TAKING ONLY FEW SAMPLES
M1<- as.matrix(scale(M_AUC))

# merge with mine
M <- Build_matrix_for_multiple_cancer_types()$Mat
#M <- Build_matrix_for_multiple_cancer_types(c(7))$Mat
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
heatmap(result)

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

result_10 <- test_all_cluster_algo(10)

result_10$fig

result_10$result

saveWidget(result_10$fig, "output/aygalic/CLUSTERING_COMPARAISON_TABLE.html", selfcontained = F, libdir = "lib")

######## ######## ######## ######## 
######## ######## ######## ######## 
######## ALIGNING CLUSTERS ######## 
######## ######## ######## ######## 
######## ######## ######## ######## 


# This function is designed to match the centroid of the 2 groups of clusters
# in a given projection space
#
# Match clusters2 to clusters1
#
# returns a new list of clusters based on clusters2, 
# but with indexes swapped to match clusters1's centroids
align_clusters_in_space <- function(clusters1, clusters2, projected_obs){
  n_clusters <- max(clusters1)
  if(n_clusters<2) return(clusters2)
  
  new_clusters <- clusters2
  for(i in 1:n_clusters){
    # find centroid of each cluster for the first dataset
    data <- projected_obs[clusters1==i,]
    centroid1 <- rapply(data, mean)
    
    # find the matching centroid
    # start with clust 1 ofc
    k = 1
    centroid2 <- rapply(projected_obs[clusters2==k,], mean)
    # compute the original distance
    d = dist(rbind(centroid1, centroid2))
    
    
    for(j in 2:n_clusters){
      centroid2_ <- rapply(projected_obs[clusters2==j,], mean)
      # compute the distance bewteen current centroids
      d_ = dist(rbind(centroid1, centroid2_))
      if(d_ < d){
        d <- d_
        centroid2 <- centroid2_
        k = j
      }
    }
    # set the new cluster k of clusters2 to the cluster i of clusters1
    new_clusters[clusters2==k] <- i
  }
  return (new_clusters)
}




# testing this function :
reduced_M1_scaled <- create_reduced_mat(M1_)


k=10
result.AUC <- kmeans(M1, centers=k)$cluster[indexes1]
# let's mess around with the clusters
new_clust <- result.AUC
# swap a few clusters
new_clust[result.AUC==1] <- 4
new_clust[result.AUC==4] <- 1

new_clust[result.AUC==8] <- 2
new_clust[result.AUC==2] <- 8

sum(!new_clust==result.AUC)

# realign ? 
new_clust <- align_clusters_in_space(result.AUC, new_clust, reduced_M1_scaled)
sum(!new_clust==result.AUC) #ok

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
v1_M1 = reduced_M1_scaled$v1
v2_M1 = reduced_M1_scaled$v2

v1_M2 = reduced_M2_scaled$v1
v2_M2 = reduced_M2_scaled$v2

names = rownames(reduced_M1_scaled)
size = length(v1_M1)

# We will use the projection axis of the RPKM dataset instead of the AUC one. 
# This shouldn't make a difference







# we prepare all the data
df1 <- data.frame(x = list(), y = list(), frame = list(), clust = list(), algo = list())
df2 <- data.frame(x = list(), y = list(), frame = list(), clust = list(), algo = list())
df3 <- data.frame(x = list(), y = list(), frame = list(), clust = list(), algo = list())
df4 <- data.frame(x = list(), y = list(), frame = list(), clust = list(), algo = list())

ALIGN = TRUE

for(i in 1:k){
  clusters1 = kmeans(M1_, centers=i)$cluster
  clusters2 = kmeans(M2_, centers=i)$cluster
  clusters3 = kmeans(M1_, centers=i)$cluster
  clusters4 = kmeans(M2_, centers=i)$cluster
  
  # make them chars
  #clusters1 <- as.character(clusters1)
  #clusters2 <- as.character(clusters2)
  #clusters3 <- as.character(clusters3)
  #clusters4 <- as.character(clusters4)
  
  if(ALIGN){
    clusters2 <- align_clusters_in_space(clusters1, clusters2, reduced_M1_scaled)
    clusters4 <- align_clusters_in_space(clusters3, clusters4, reduced_M2_scaled)
  }
  
  df1 <- rbind(df1, data.frame(x = v1_M1, y = v2_M1, frame = rep(i, size), name = names,
                               clust = clusters1, 
                               algo = rep(all_algos[[1]], size)))
  
  df2 <- rbind(df2, data.frame(x = v1_M1, y = v2_M1, frame = rep(i, size), name = names,
                               clust = clusters2, 
                               algo = rep(all_algos[[1]], size)))
  
  df3 <- rbind(df3, data.frame(x = v1_M2, y = v2_M2, frame = rep(i, size), name = names,
                               clust = clusters3, 
                               algo = rep(all_algos[[1]], size)))
  
  df4 <- rbind(df4, data.frame(x = v1_M2, y = v2_M2, frame = rep(i, size), name = names,
                               clust = clusters4, 
                               algo = rep(all_algos[[1]], size)))
}


for(j in 2:n_algo){
  for(i in 1:k){
    clusters1 = hcut(M1_, i, hc_method = algos[[j-1]][2], hc_metric= algos[[j-1]][1])$cluster
    clusters2 = hcut(M2_, i, hc_method = algos[[j-1]][2], hc_metric= algos[[j-1]][1])$cluster
    clusters3 = hcut(M1_, i, hc_method = algos[[j-1]][2], hc_metric= algos[[j-1]][1])$cluster
    clusters4 = hcut(M2_, i, hc_method = algos[[j-1]][2], hc_metric= algos[[j-1]][1])$cluster
    
    # make them chars
    #clusters1 <- as.character(clusters1)
    #clusters2 <- as.character(clusters2)
    #clusters3 <- as.character(clusters3)
    #clusters4 <- as.character(clusters4)
    
    if(ALIGN){
      clusters2 <- align_clusters_in_space(clusters1, clusters2, reduced_M1_scaled)
      clusters4 <- align_clusters_in_space(clusters3, clusters4, reduced_M2_scaled)
    }
    
    df1 <- rbind(df1, data.frame(x = v1_M1, y = v2_M1, frame = rep(i, size), name = names,
                                 clust = clusters1, 
                                 algo = rep(all_algos[[j]], size)))
    
    df2 <- rbind(df2, data.frame(x = v1_M1, y = v2_M1, frame = rep(i, size), name = names,
                                 clust = clusters2,
                                 algo = rep(all_algos[[j]], size)))
    
    df3 <- rbind(df3, data.frame(x = v1_M2, y = v2_M2, frame = rep(i, size), name = names,
                                 clust = clusters3, 
                                 algo = rep(all_algos[[j]], size)))
    
    df4 <- rbind(df4, data.frame(x = v1_M2, y = v2_M2, frame = rep(i, size), name = names,
                                 clust = clusters4,
                                 algo = rep(all_algos[[j]], size)))
    
    
  }
}




# we build the plots
fig1 <- plot_ly()
fig2 <- plot_ly()
fig3 <- plot_ly()
fig4 <- plot_ly()



for(alg in all_algos){
  # optimization
  VISIBLE = ifelse(alg=="kmeans", T, F)
  df1_ <- df1[df1$algo == alg,]
  df2_ <- df2[df2$algo == alg,]
  df3_ <- df3[df3$algo == alg,]
  df4_ <- df4[df4$algo == alg,]
  
  fig1 <- fig1 %>% add_markers(
    data = df1_,
    x = ~x, y = ~y,
    text = paste(df1_$name, alg),
    hoverinfo = "text",
    frame = ~frame,
    color = ~clust,
    marker =list(colorscale = 'Jet'),
    showlegend = F, 
    visible = VISIBLE
  )
  fig2 <- fig2 %>% add_markers(
    data = df2_,
    x = ~x, y = ~y,
    text = paste(df2_$name, alg),
    hoverinfo = "text",
    frame = ~frame,
    color = ~clust,
    marker =list(colorscale = 'Jet'),
    showlegend = F,
    visible = VISIBLE
  )
  fig3 <- fig3 %>% add_markers(
    data = df3_,
    x = ~x, y = ~y,
    text = paste(df3_$name, alg),
    hoverinfo = "text",
    frame = ~frame,
    color = ~clust,
    marker =list(colorscale = 'Jet'),
    showlegend = F, 
    visible = VISIBLE
  )
  fig4 <- fig4 %>% add_markers(
    data = df4_,
    x = ~x, y = ~y,
    text = paste(df4_$name, alg),
    hoverinfo = "text",
    frame = ~frame,
    color = ~clust,
    marker =list(colorscale = 'Jet'),
    showlegend = F,
    visible = VISIBLE
  )
}

# add algorithm selection
fig <- subplot(fig1, fig2, fig3, fig4, nrows = 2, which_layout = 0) 

# We are showing 4 traces at once : you can see this "matrix" in the following way:
# 4 square matrices with a diagonal "TRUE" separated by a column of "FALSE"
# extra column of "FALSE" is added at the end
# The columns of false is here to make everything coherent
METHOD = "restlye"

BTN2 = list(  
  list(method = METHOD, args = list(list(visible = c(T, F, F, F, F, F, F, F, F, F, F, T, F, F, F, F, F, F, F, F, F, F, T, F, F, F, F, F, F, F, F, F, F, T, F, F, F, F, F, F, F, F, F, F))), label = "kmeans"),
  list(method = METHOD, args = list(list(visible = c(F, T, F, F, F, F, F, F, F, F, F, F, T, F, F, F, F, F, F, F, F, F, F, T, F, F, F, F, F, F, F, F, F, F, T, F, F, F, F, F, F, F, F, F))), label = "euclidian single"),       
  list(method = METHOD, args = list(list(visible = c(F, F, T, F, F, F, F, F, F, F, F, F, F, T, F, F, F, F, F, F, F, F, F, F, T, F, F, F, F, F, F, F, F, F, F, T, F, F, F, F, F, F, F, F))), label = "euclidian average"),        
  list(method = METHOD, args = list(list(visible = c(F, F, F, T, F, F, F, F, F, F, F, F, F, F, T, F, F, F, F, F, F, F, F, F, F, T, F, F, F, F, F, F, F, F, F, F, T, F, F, F, F, F, F, F))), label = "euclidian complete"),        
  list(method = METHOD, args = list(list(visible = c(F, F, F, F, T, F, F, F, F, F, F, F, F, F, F, T, F, F, F, F, F, F, F, F, F, F, T, F, F, F, F, F, F, F, F, F, F, T, F, F, F, F, F, F))), label = "manhattan single"),        
  list(method = METHOD, args = list(list(visible = c(F, F, F, F, F, T, F, F, F, F, F, F, F, F, F, F, T, F, F, F, F, F, F, F, F, F, F, T, F, F, F, F, F, F, F, F, F, F, T, F, F, F, F, F))), label = "manhattan average"),        
  list(method = METHOD, args = list(list(visible = c(F, F, F, F, F, F, T, F, F, F, F, F, F, F, F, F, F, T, F, F, F, F, F, F, F, F, F, F, T, F, F, F, F, F, F, F, F, F, F, T, F, F, F, F))), label = "manhattan complete"),        
  list(method = METHOD, args = list(list(visible = c(F, F, F, F, F, F, F, T, F, F, F, F, F, F, F, F, F, F, T, F, F, F, F, F, F, F, F, F, F, T, F, F, F, F, F, F, F, F, F, F, T, F, F, F))), label = "canberra single"),        
  list(method = METHOD, args = list(list(visible = c(F, F, F, F, F, F, F, F, T, F, F, F, F, F, F, F, F, F, F, T, F, F, F, F, F, F, F, F, F, F, T, F, F, F, F, F, F, F, F, F, F, T, F, F))), label = "canberra average"),        
  list(method = METHOD, args = list(list(visible = c(F, F, F, F, F, F, F, F, F, T, F, F, F, F, F, F, F, F, F, F, T, F, F, F, F, F, F, F, F, F, F, T, F, F, F, F, F, F, F, F, F, F, T, F))), label = "canberra complete")
)



fig <- fig %>% animation_opts(1000, easing = "elastic", redraw = F) 
fig <- fig %>% layout(title = "Comparing clustering Algo between the 2 datasets",
                      updatemenus = list(list(y = 0.8, active = 0, type= 'buttons', buttons = BTN2)))
fig <- fig %>% layout(plot_bgcolor='#ddd',
                      xaxis  = list(zeroline = F, showgrid = F, showticklabels=FALSE),
                      xaxis2 = list(zeroline = F, showgrid = F, showticklabels=FALSE),
                      xaxis3 = list(zeroline = F, showgrid = F, showticklabels=FALSE),
                      xaxis4 = list(zeroline = F, showgrid = F, showticklabels=FALSE),
                      yaxis  = list(zeroline = F, showgrid = F, showticklabels=FALSE),
                      yaxis2 = list(zeroline = F, showgrid = F, showticklabels=FALSE),
                      yaxis3 = list(zeroline = F, showgrid = F, showticklabels=FALSE),
                      yaxis4 = list(zeroline = F, showgrid = F, showticklabels=FALSE),
                      showlegend = F)
                        
# doesn't work
#fig %>% hide_colorbar()

# add javascript because plotly is an undocumented glitchy sesspool
# fig %>% onRender("")

fig



saveWidget(fig, "output/aygalic/CLUSTERING_COMPARAISON.html", selfcontained = F, libdir = "lib")

#############################################
#############################################
###### ONE LAST PLOT I WANT TO DO HERE ######
#############################################
#############################################
# I want to do just like the first plot where
# i compare all the clustering algo but this
# time i compare them to the primary tissue
# 
# So we will endup with two tables, one from
# the RPKM dataset, another one from the AUC
# dataset.
#############################################
#############################################
#############################################

primitive_types <- str_split(names, '_', n = 2, simplify = TRUE)[,2]

k=10
result.AUC <- kmeans(M1, centers=k)$cluster[indexes1]
result.rpkm <- kmeans(M2, centers=k)$cluster[indexes2]

rownames(primitive_types) <- rownames(result.AUC) 


result1 <- table(result.AUC, primitive_types)
result2 <- table(result.rpkm, primitive_types)

result1
fig1 <- plot_ly(z = result1, type = "heatmap")
fig2 <- plot_ly(z = result2, type = "heatmap")
fig <- subplot(fig1, fig2)
fig

saveWidget(fig, "output/aygalic/clutsers_primitive_comparaison.html", selfcontained = F, libdir = "lib")



# i lied i'm not done

# second try at this
# ONLY WORKS IN R^2 FOR NOW
align_clusters_in_space2 <- function(clusters1, clusters2, projected_obs){
  n_clusters <- max(clusters1)
  if(n_clusters<2) return(clusters2)
  
  centroids1 <- data.frame(x=list(), y=list())
  centroids2 <- data.frame(x=list(), y=list())
  

  
  for(i in 1:n_clusters){
    data <- projected_obs[clusters1==i,]
    mu = rapply(data, mean)
    centroids1[i, 1] <- mu[1]
    centroids1[i, 2] <- mu[2]
  }
  for(i in 1:n_clusters){
    data <- projected_obs[clusters2==i,]
    mu = rapply(data, mean)
    centroids2[i, 1] <- mu[1]
    centroids2[i, 2] <- mu[2]  
  }
  colnames(centroids1) <- c("x","y")
  colnames(centroids2) <- c("x","y")
  
  # we want to "align" the centroids
  # Let's first compute the distance 
  f_dist <- function(a, b){
    return(sqrt(((a$x - b$x)^2)+(a$y - b$y)^2))
  } 
    
  dist.e = data.frame()
  for(i in 1:n_clusters){
    for(j in 1:n_clusters){
      dist.e[i,j] = f_dist(centroids1[i,], centroids2[j,])
    }
  }
  
  # let's now do the ranking and ordering
  order = rapply(dist.e, which.min)
  
  new_clust <- as.factor(clusters1)
  levels(new_clust) <- order
  
  
  
  return(list(centroids1_before = centroids1, 
              centroids2_before = centroids2, 
              dist = dist.e, 
              order = order, 
              #new_centroid = new_centroid,
              new_clusters = new_clust ))
}


k=3
c1 <- kmeans(M1, centers=k)$cluster[indexes1]
c2 <- kmeans(M2, centers=k)$cluster[indexes2]






test <- align_clusters_in_space2(c1, c2, reduced_M2_scaled)
test

as.factor(c1)==test$new_clusters


plot(reduced_M2_scaled, col = c1, pch = 20, lwd = 0)

points(test$centroids1_before, col = "blue", pch=3, lwd = 5)
points(test$centroids2_before, col = "yellow", pch=3, lwd = 5)

points(test$centroids1_before[1,], col = "green", pch=3, lwd = 5)
points(test$centroids2_before[1,], col = "green", pch=3, lwd = 5)





plot(test)

c2 <- align_clusters_in_space(c1, c2, reduced_M2_scaled)

resulttable <- table(c1, c2)
resulttable
plot_ly(z = resulttable, type = "heatmap")


table(result.rpkm, result.AUC)

plot_ly(z = table(result.rpkm, result.AUC), type = "heatmap")
