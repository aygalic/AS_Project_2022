library(ggplot2)
library(plotly)  # interactive plots 
library(mvtnorm)
library(rgl)
library(car)

setwd("~/OneDrive/polimi/COURSES/S8/APPLIED_STATS/AS_Project_2022")
# all the useful functions for pca and projections
source("src/utilities.R")



#import data
data_patient_= read.delim(file.path("Dataset", "data_clinical_patient.txt"), header = TRUE, comment.char = '#')
data_sample_= read.delim(file.path("Dataset",'data_clinical_sample.txt'), header = TRUE, comment.char = '#')
# use processed data
original_data_mrna_ = read.delim(file.path("Dataset", "1_rpkm.txt"), header = TRUE, comment.char = '#', nrows=5000)
# we might as well use all the data
original_data_mrna_ = read.delim(file.path("Dataset", "1_rpkm.txt"), header = TRUE, comment.char = '#')





# Get a table of all cancer types and occurrences
cancer_types <- as.data.frame(table(data_sample_$"CANCER_TYPE_DETAILED"), stringsAsFactors = FALSE)
cancer_types <- cancer_types[order(cancer_types$Freq, decreasing = TRUE),]
names(cancer_types)<-c("Factor", "Freq")




# creating a matrix will all cancer types
#
# We could also create a matrix with specified cancer types only
cancer_type_selection <- c(1:(length(cancer_types$Factor)-1))
M = Build_matrix_for_multiple_cancer_types(cancer_type_selection)
tags = M$Tags
M <- M$Mat
M_scaled <- as.matrix(scale(M))



M_ <- M_scaled

# We want to cluster cells together so we use the transpose
M_ <-t(M)


M_.e <- dist(M_, method='euclidean')
M_.m <- dist(M_, method='manhattan')
M_.c <- dist(M_, method='canberra')



# making every possible tree
M_.es <- hclust(M_.e, method='single')
M_.ea <- hclust(M_.e, method='average')
M_.ec <- hclust(M_.e, method='complete')

M_.ms <- hclust(M_.m, method='single')
M_.ma <- hclust(M_.m, method='average')
M_.mc <- hclust(M_.m, method='complete')

M_.cs <- hclust(M_.c, method='single')
M_.ca <- hclust(M_.c, method='average')
M_.cc <- hclust(M_.c, method='complete')





# We generate vectors of labels through the command cutree()
help(cutree)

# Fix k=2 clusters:
cluster.ec <- cutree(M_.ec, k=2) # euclidean-complete:
cluster.ec


library(factoextra) # clustering visualization

library(ggpubr) # ggarrange

# Making a silhouette plot with all metods
p1 <- fviz_nbclust(M_, FUN = hcut, method = "silhouette", hc_method ="complete", hc_metric="euclidian")
p2 <- fviz_nbclust(M_, FUN = hcut, method = "silhouette", hc_method ="complete", hc_metric="manhattan")
p3 <- fviz_nbclust(M_, FUN = hcut, method = "silhouette", hc_method ="complete", hc_metric="canberra")

p4 <- fviz_nbclust(M_, FUN = hcut, method = "silhouette", hc_method ="average", hc_metric="euclidian")
p5 <- fviz_nbclust(M_, FUN = hcut, method = "silhouette", hc_method ="average", hc_metric="manhattan")
p6 <- fviz_nbclust(M_, FUN = hcut, method = "silhouette", hc_method ="average", hc_metric="canberra")

p7 <- fviz_nbclust(M_, FUN = hcut, method = "silhouette", hc_method ="single", hc_metric="euclidian")
p8 <- fviz_nbclust(M_, FUN = hcut, method = "silhouette", hc_method ="single", hc_metric="manhattan")
p9 <- fviz_nbclust(M_, FUN = hcut, method = "silhouette", hc_method ="single", hc_metric="canberra")

silhouette_plot_all_method <- ggarrange(p1, p2, p3, p4, p5, p6, p7, p8 , p9,
                                        ncol = 3, nrow = 3)
x11()
silhouette_plot_all_method


# Making a ELBOW plot with all metods
e1 <- fviz_nbclust(M_, FUN = hcut, method = "wss", hc_method ="complete", hc_metric="euclidian")
e2 <- fviz_nbclust(M_, FUN = hcut, method = "wss", hc_method ="complete", hc_metric="manhattan")
e3 <- fviz_nbclust(M_, FUN = hcut, method = "wss", hc_method ="complete", hc_metric="canberra")

e4 <- fviz_nbclust(M_, FUN = hcut, method = "wss", hc_method ="average", hc_metric="euclidian")
e5 <- fviz_nbclust(M_, FUN = hcut, method = "wss", hc_method ="average", hc_metric="manhattan")
e6 <- fviz_nbclust(M_, FUN = hcut, method = "wss", hc_method ="average", hc_metric="canberra")

e7 <- fviz_nbclust(M_, FUN = hcut, method = "wss", hc_method ="single", hc_metric="euclidian")
e8 <- fviz_nbclust(M_, FUN = hcut, method = "wss", hc_method ="single", hc_metric="manhattan")
e9 <- fviz_nbclust(M_, FUN = hcut, method = "wss", hc_method ="single", hc_metric="canberra")

elbow_plot_all_method <- ggarrange(e1, e2, e3, e4, e5, e6, e7, e8 , e9,
                                        ncol = 3, nrow = 3)
x11()
elbow_plot_all_method









# WHAT HAPPENS IF YOU PLOT THOSES ON THE REDUCED SPACE PROVIDED BY PCA ?
# create a projected space 
reduced_M_scaled = create_reduced_mat(M_, 2)

 # big brain graph

# optimal number of cluster would be between 4-7
k=7
cluster.es <- cutree(M_.es, k=k)
cluster.ea <- cutree(M_.ea, k=k)
cluster.ec <- cutree(M_.ec, k=k)

cluster.ms <- cutree(M_.ms, k=k)
cluster.ma <- cutree(M_.ma, k=k)
cluster.mc <- cutree(M_.mc, k=k)

cluster.cs <- cutree(M_.cs, k=k)
cluster.ca <- cutree(M_.ca, k=k)
cluster.cc <- cutree(M_.cc, k=k)

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




# Big brain plot

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
v1 = reduced_M_scaled$v1
v2 = reduced_M_scaled$v2
names = rownames(reduced_M_scaled)
size = length(v1)
k = 10


# FIRST WE BUILD THE FIRST TRACE 
# We une this one with kmeans since it's the only different algo

df <- data.frame(x = list(), y = list(), frame = list(), clust = list())
for(i in 1:k){
  df <- rbind(df, data.frame(x = v1, y = v2, frame = rep(i, size), name = names,
                             clust = kmeans(M_, centers=i)$cluster))
}

fig <- df %>% plot_ly()

fig <- fig %>% add_markers(
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
  df <- data.frame( x = list(), y = list(), frame = list(), clust = list())
  for(i in 1:k){
    df <- rbind(df, data.frame(x = v1, y = v2, frame = rep(i ,size), name = names,
                               clust = hcut(M_, i, hc_method = algo[2], hc_metric= algo[1])$cluster))
  }
  
  fig <- fig %>% add_markers(
    data = df, x = ~x, y = ~y,
    text = ~name,
    hoverinfo = "text",
    frame = ~frame,
    color = ~clust,
    marker =list(colorscale = 'Viridis'),
    showlegend = F,
    visible = F
  )
}



# add algorithm selection
fig <- fig %>% layout(
  title = "Comparing clustering Algo",
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


saveWidget(fig, "output/aygalic/Clust_cells.html", selfcontained = F, libdir = "lib")



# same with breast cancer

#

cancer_type_selection <- c(7)
M = Build_matrix_for_multiple_cancer_types(cancer_type_selection)
tags = M$Tags
M <- M$Mat
M_ <- as.matrix(scale(M))
M_ <-t(M)
reduced_M_scaled = create_reduced_mat(M_, 2)

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
v1 = reduced_M_scaled$v1
v2 = reduced_M_scaled$v2
names = rownames(reduced_M_scaled)
size = length(v1)
k = 10


# FIRST WE BUILD THE FIRST TRACE 
# We une this one with kmeans since it's the only different algo

df <- data.frame(x = list(), y = list(), frame = list(), clust = list())
for(i in 1:k){
  df <- rbind(df, data.frame(x = v1, y = v2, frame = rep(i, size), name = names,
                             clust = kmeans(M_, centers=i)$cluster))
}

fig <- df %>% plot_ly()

fig <- fig %>% add_markers(
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
  df <- data.frame( x = list(), y = list(), frame = list(), clust = list())
  for(i in 1:k){
    df <- rbind(df, data.frame(x = v1, y = v2, frame = rep(i ,size), name = names,
                               clust = hcut(M_, i, hc_method = algo[2], hc_metric= algo[1])$cluster))
  }
  
  fig <- fig %>% add_markers(
    data = df, x = ~x, y = ~y,
    text = ~name,
    hoverinfo = "text",
    frame = ~frame,
    color = ~clust,
    marker =list(colorscale = 'Viridis'),
    showlegend = F,
    visible = F
  )
}



# add algorithm selection
fig <- fig %>% layout(
  title = "Comparing clustering Algo for BREAST CANCER ONLY",
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


saveWidget(fig, "output/aygalic/Clust_cells_BREAST.html", selfcontained = F, libdir = "lib")


