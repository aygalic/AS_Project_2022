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





# this is a plot of how similar are the classifications between both datasets 
# depending on the number of cluster

result_10 <- test_all_cluster_algo(10)

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

# AUC matrix
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


# AUC projections
v1_M2 = reduced_M2_scaled$v1
v2_M2 = reduced_M2_scaled$v2

names = rownames(reduced_M2_scaled)
size = length(v1_M2)



# we prepare all the data
df3 <- data.frame(x = list(), y = list(), frame = list(), clust = list(), algo = list())
df4 <- data.frame(x = list(), y = list(), frame = list(), clust = list(), algo = list())

ALIGN = TRUE

for(i in 1:k){
  clusters3 = kmeans(M1_, centers=i)$cluster
  clusters4 = kmeans(M2_, centers=i)$cluster
  
  # make them chars
  #clusters1 <- as.character(clusters1)
  #clusters2 <- as.character(clusters2)
  #clusters3 <- as.character(clusters3)
  #clusters4 <- as.character(clusters4)
  
  if(ALIGN){
    clusters4 <- align_clusters_in_space(clusters3, clusters4, reduced_M2_scaled)
  }

  
  df3 <- rbind(df3, data.frame(x = v1_M2, y = v2_M2, frame = rep(i, size), name = names,
                               clust = clusters3, 
                               algo = rep(all_algos[[1]], size)))
  
  df4 <- rbind(df4, data.frame(x = v1_M2, y = v2_M2, frame = rep(i, size), name = names,
                               clust = clusters4, 
                               algo = rep(all_algos[[1]], size)))
}


for(j in 2:n_algo){
  for(i in 1:k){
    clusters3 = hcut(M1_, i, hc_method = algos[[j-1]][2], hc_metric= algos[[j-1]][1])$cluster
    clusters4 = hcut(M2_, i, hc_method = algos[[j-1]][2], hc_metric= algos[[j-1]][1])$cluster
    
    # make them chars
    #clusters1 <- as.character(clusters1)
    #clusters2 <- as.character(clusters2)
    #clusters3 <- as.character(clusters3)
    #clusters4 <- as.character(clusters4)
    
    if(ALIGN){
      clusters4 <- align_clusters_in_space(clusters3, clusters4, reduced_M2_scaled)
    }
    
    
    df3 <- rbind(df3, data.frame(x = v1_M2, y = v2_M2, frame = rep(i, size), name = names,
                                 clust = clusters3, 
                                 algo = rep(all_algos[[j]], size)))
    
    df4 <- rbind(df4, data.frame(x = v1_M2, y = v2_M2, frame = rep(i, size), name = names,
                                 clust = clusters4,
                                 algo = rep(all_algos[[j]], size)))
    
    
  }
}




# we build the plots
fig3 <- plot_ly()
fig4 <- plot_ly()



for(alg in all_algos){
  # optimization
  VISIBLE = ifelse(alg=="kmeans", T, F)
  df3_ <- df3[df3$algo == alg,]
  df4_ <- df4[df4$algo == alg,]
  

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
fig <- subplot(fig4, fig3, nrows = 1, which_layout = 0) 

# We are showing 4 traces at once : you can see this "matrix" in the following way:
# 4 square matrices with a diagonal "TRUE" separated by a column of "FALSE"
# extra column of "FALSE" is added at the end
# The columns of false is here to make everything coherent
METHOD = "restlye"

BTN2 = list(  
  list(method = METHOD, args = list(list(visible = c(T, F, F, F, F, F, F, F, F, F, F, T, F, F, F, F, F, F, F, F, F, F))), label = "kmeans"),
  list(method = METHOD, args = list(list(visible = c(F, T, F, F, F, F, F, F, F, F, F, F, T, F, F, F, F, F, F, F, F, F))), label = "euclidian single"),       
  list(method = METHOD, args = list(list(visible = c(F, F, T, F, F, F, F, F, F, F, F, F, F, T, F, F, F, F, F, F, F, F))), label = "euclidian average"),        
  list(method = METHOD, args = list(list(visible = c(F, F, F, T, F, F, F, F, F, F, F, F, F, F, T, F, F, F, F, F, F, F))), label = "euclidian complete"),        
  list(method = METHOD, args = list(list(visible = c(F, F, F, F, T, F, F, F, F, F, F, F, F, F, F, T, F, F, F, F, F, F))), label = "manhattan single"),        
  list(method = METHOD, args = list(list(visible = c(F, F, F, F, F, T, F, F, F, F, F, F, F, F, F, F, T, F, F, F, F, F))), label = "manhattan average"),        
  list(method = METHOD, args = list(list(visible = c(F, F, F, F, F, F, T, F, F, F, F, F, F, F, F, F, F, T, F, F, F, F))), label = "manhattan complete"),        
  list(method = METHOD, args = list(list(visible = c(F, F, F, F, F, F, F, T, F, F, F, F, F, F, F, F, F, F, T, F, F, F))), label = "canberra single"),        
  list(method = METHOD, args = list(list(visible = c(F, F, F, F, F, F, F, F, T, F, F, F, F, F, F, F, F, F, F, T, F, F))), label = "canberra average"),        
  list(method = METHOD, args = list(list(visible = c(F, F, F, F, F, F, F, F, F, T, F, F, F, F, F, F, F, F, F, F, T, F))), label = "canberra complete")
)



fig <- fig %>% animation_opts(1000, easing = "elastic", redraw = F) 
fig <- fig %>% layout(title = "Comparing clustering Algo between the 2 datasets",
                      updatemenus = list(list(y = 0.8, active = 0, type= 'buttons', buttons = BTN2)))
fig <- fig %>% layout(plot_bgcolor='#ddd',
                      xaxis  = list(zeroline = F, showgrid = F, showticklabels=FALSE),
                      xaxis2 = list(zeroline = F, showgrid = F, showticklabels=FALSE),
                      yaxis  = list(zeroline = F, showgrid = F, showticklabels=FALSE),
                      yaxis2 = list(zeroline = F, showgrid = F, showticklabels=FALSE),
                      showlegend = F)
                        
# doesn't work
#fig %>% hide_colorbar()

# add javascript because plotly is an undocumented glitchy sesspool
# fig %>% onRender("")

fig




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

# we have to normalize based on cancer count
# but first we order them
reorder <- function(mat, rescale = T){
  sum = colSums(mat)
  if(!rescale)
    return ( mat[, order(sum, decreasing = T)])
  mat <- apply(mat, 2, scale)
  return ( mat[, order(sum, decreasing = T)])
}


primitive_types <- str_split(names, '_', n = 2, simplify = TRUE)[,2]

k=10
result.AUC <- kmeans(M1, centers=k)$cluster[indexes1]
result.rpkm <- kmeans(M2, centers=k)$cluster[indexes2]

rownames(primitive_types) <- rownames(result.AUC) 





result1 <- table(result.AUC, primitive_types)
result2 <- table(result.rpkm, primitive_types)

result1 <- reorder(result1, F)
result2 <- reorder(result2, F)

as.matrix(result1)
as.data.frame(as.matrix(result1))
(result1)

fig1 <- plot_ly(z = result1, type = "heatmap", x = colnames(result1),  y = rownames(result1))          
fig2 <- plot_ly(z = result2, type = "heatmap", x = colnames(result2),  y = rownames(result2))
fig <- subplot(fig1, fig2)
fig



# we prepare all the data
df1 <- data.frame()
df2 <- data.frame()

ALIGN = TRUE

for(i in 2:k){
  # should I use M1_ and M2_ instead ?
  result.AUC <- kmeans(M1, centers=i)$cluster[indexes1]
  result.rpkm <- kmeans(M2, centers=i)$cluster[indexes2]
  result1 <- table(result.AUC, primitive_types)
  result2 <- table(result.rpkm, primitive_types)
  
  result1 <- as.data.frame(reorder(result1, F))
  result2 <- as.data.frame(reorder(result2, F))
  
  result1$algo = rep(all_algos[[1]], i)
  result2$algo = rep(all_algos[[1]], i)
  
  result1$frame = rep(i, i)
  result2$frame = rep(i, i)
  print(i)
  print(result1)
  df1 <- rbind(df1, result1)
  df2 <- rbind(df2, result2)
}


for(j in 2:n_algo){
  for(i in 1:k){
    result.AUC = hcut(M1_, i, hc_method = algos[[j-1]][2], hc_metric= algos[[j-1]][1])$cluster
    result.rpkm = hcut(M2_, i, hc_method = algos[[j-1]][2], hc_metric= algos[[j-1]][1])$cluster
    
    result1 <- table(result.AUC, primitive_types)
    result2 <- table(result.rpkm, primitive_types)
    
    result1 <- as.data.frame(reorder(result1, F))
    result2 <- as.data.frame(reorder(result2, F))
    
    result1$algo = rep(all_algos[[1]], i)
    result2$algo = rep(all_algos[[1]], i)
    
    result1$frame = rep(i, i)
    result2$frame = rep(i, i)
    
    df1 <- rbind(df1, result1)
    df2 <- rbind(df2, result2)
  }
}

# we build the plots
fig1 <- plot_ly()
fig2 <- plot_ly()

for(alg in all_algos){
  # optimization
  VISIBLE = ifelse(alg=="kmeans", T, F)
  df1_ <- df1[df1$algo == alg,]
  df2_ <- df1[df1$algo == alg,]
  
  fig1 <- plot_ly(z = result1, type = "heatmap", x = colnames(result1),  y = rownames(result1), frame = result1$frame)          
  fig2 <- plot_ly(z = result2, type = "heatmap", x = colnames(result2),  y = rownames(result2), frame = result2$frame)
  
  #fig1 <- fig1 %>% add_markers(
  #  data = df3_,
  #  x = ~x, y = ~y,
  #  text = paste(df3_$name, alg),
  #  hoverinfo = "text",
  #  frame = ~frame,
  #  color = ~clust,
  #  marker =list(colorscale = 'Jet'),
  #  showlegend = F, 
  #  visible = VISIBLE
  #)
  #fig2 <- fig2 %>% add_markers(
  #  data = df4_,
  #  x = ~x, y = ~y,
  #  text = paste(df4_$name, alg),
  #  hoverinfo = "text",
  #  frame = ~frame,
  #  color = ~clust,
  #  marker =list(colorscale = 'Jet'),
  #  showlegend = F,
  #  visible = VISIBLE
  #)
}

# add algorithm selection
fig <- subplot(fig1, fig2, nrows = 1, which_layout = 0) 



















