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
original_data_mrna_ = read.delim(file.path("Dataset", "1_rpkm.txt"), header = TRUE, comment.char = '#')

# Get a table of all cancer types and occurrences
cancer_types <- as.data.frame(table(data_sample_$"CANCER_TYPE_DETAILED"), stringsAsFactors = FALSE)
cancer_types <- cancer_types[order(cancer_types$Freq, decreasing = TRUE),]
names(cancer_types)<-c("Factor", "Freq")


# creating a matrix with all cancer types
#
# We could also create a matrix with specified cancer types only
cancer_type_selection <- c(1:(length(cancer_types$Factor)-1))
M = Build_matrix_for_multiple_cancer_types(cancer_type_selection)
tags = M$Tags
M <- M$Mat
M_scaled <- as.matrix(scale(M))




### in automatic, command kmeans()
help(kmeans)

result.k <- kmeans(M_scaled, centers=2) # Centers: fixed number of clusters



# WHAT HAPPENS IF YOU PLOT THOSES ON THE REDUCED SPACE PROVIDED BY PCA ?
# create a projected space 
reduced_M_scaled = create_reduced_mat(M_scaled, 2)


p1 <- plot_ly(data = data.frame(reduced_M_scaled), x = ~v1, y = ~v2,
        text = rownames(reduced_M_scaled), type = "scatter", 
        color = factor(result.k$cluster)) %>% 
  layout(margin = c(10,10,10,10,0))
p1 




# experimenting with the number of clusters
result.k <- kmeans(M_scaled, centers=3) # Centers: fixed number of clusters
p2 <- plot_ly(data = data.frame(reduced_M_scaled), x = ~v1, y = ~v2,
        text = rownames(reduced_M_scaled), type = "scatter", 
        color = factor(result.k$cluster)) %>% 
  layout(margin = c(10,10,10,10,0))
p2

# experimenting with the number of clusters
result.k <- kmeans(M_scaled, centers=5) # Centers: fixed number of clusters
plot_ly(data = data.frame(reduced_M_scaled), x = ~v1, y = ~v2,
        text = rownames(reduced_M_scaled), type = "scatter", 
        color = factor(result.k$cluster)) %>% 
  layout(margin = c(10,10,10,10,0))

plot_ly(data = data.frame(reduced_M_scaled), x = ~v1, y = ~v2,
        mode = "markers", 
        marker = list(color = factor(result.k$cluster)),
        text = rownames(reduced_M_scaled), type = "scatter") %>% 
  layout(margin = c(10,10,10,10,0))






# BIGG BRAIN PLOT



j = 2
k = 10
aval <- list()
for(step in j:k){
  aval[[step]] <-list(visible = FALSE,
                      x = reduced_M_scaled$v1,
                      y = reduced_M_scaled$v2,
                      name = step,
                      col = kmeans(M_scaled, centers=step)$cluster)
}
aval[3][[1]]$visible = TRUE

# create steps and plot all traces
steps <- list()
fig <- plot_ly()
for (i in j:k) {
  fig <- add_trace(fig,x=aval[i][[1]]$x,  y=aval[i][[1]]$y, visible = aval[i][[1]]$visible,
                   type = 'scatter', mode = 'markers', marker=list(color=aval[i][[1]]$col),
                   name = aval[i][[1]]$name)
  
  step <- list(args = list('visible', rep(FALSE, length(aval))),
               method = 'restyle')
  step$args[[2]][i] = TRUE  
  steps[[i]] = step 
}  

# add slider control to plot
fig <- fig %>%
  layout(sliders = list(list(active = 3,
                             currentvalue = list(prefix = "Number of clust: "),
                             steps = steps)))

fig
