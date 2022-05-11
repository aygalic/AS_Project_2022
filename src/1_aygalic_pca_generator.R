library(ggplot2)
library(plotly)  # interactive plots 

setwd("~/OneDrive/polimi/COURSES/S8/APPLIED_STATS/AS_Project_2022")

# all the useful functions for pca and projections
source("src/utilities.R")

#import data
#
# This override every data imported in utilities.R
data_patient_= read.delim(file.path("Dataset", "data_clinical_patient.txt"), header = TRUE, comment.char = '#')
data_sample_= read.delim(file.path("Dataset",'data_clinical_sample.txt'), header = TRUE, comment.char = '#')
original_data_mrna_ = read.delim(file.path("Dataset", "1_rpkm.txt"), header = TRUE, comment.char = '#', nrows=5000)

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





#############################
############ PCA ############
#############################

# We reduce the space of cell line and project genes (to see if there are any gene of interest ?)
reduced_M_scaled = create_reduced_mat(M_scaled, 4)

p1 <- plot_ly(data = data.frame(reduced_M_scaled), x = ~v1, y = ~v2,
             text = rownames(reduced_M_scaled), type = "scatter") %>% 
  layout(margin = c(10,10,10,10,0))

p1


# saving plot 
saveWidget(p1, "output/aygalic/PCA_base.html", selfcontained = F, libdir = "lib")


# adding color based on overall gene expression (before scaling)

gene_exp = apply(M, 1, sum)

p2 <- plot_ly(data = data.frame(reduced_M_scaled), x = ~v1, y = ~v2,
        text = rownames(reduced_M_scaled), type = "scatter",
        color = gene_exp) %>% 
  layout(margin = c(10,10,10,10,0))
p2

# same but log scaled ?

gene_exp = apply(M, 1, sum)
gene_exp_scaled = log(gene_exp)

p3 <- plot_ly(data = data.frame(reduced_M_scaled), x = ~v1, y = ~v2,
        text = rownames(reduced_M_scaled), type = "scatter",
        color = gene_exp_scaled) %>% 
  layout(margin = c(10,10,10,10,0))
p3
saveWidget(p3, "output/aygalic/PCA_log_scaled_colors.html", selfcontained = F, libdir = "lib")



# We try a 3D plot 
plot_ly(data = data.frame(reduced_M_scaled), x = ~v1, y = ~v2, z = ~v3,
        text = rownames(reduced_M_scaled),
        color = gene_exp_scaled) %>% 
  layout(margin = c(10,10,10,10,0))

# it's pretty bad compared to the 2D plot which gives a way better diplay of  
# the observations + it takes forever to render


# WE WANT TO INVESTIGATE THE 2 AXIS WE FOUND : 


# investigate axis 3
fig1 <- plot_ly(data = data.frame(reduced_M_scaled), x = ~v1, y = ~v3,
                text = rownames(reduced_M_scaled), type = "scatter", name = 'v1/v3')

fig2 <- plot_ly(data = data.frame(reduced_M_scaled), x = ~v2, y = ~v3,
                text = rownames(reduced_M_scaled), type = "scatter", name = 'v2/v3')

p4 <- subplot(fig1, fig2, nrows = 1) 
p4
saveWidget(p4, "output/aygalic/PCA_axis_3.html", selfcontained = F, libdir = "lib")


# investigate axis 4
fig1 <- plot_ly(data = data.frame(reduced_M_scaled), x = ~v1, y = ~v4,
                text = rownames(reduced_M_scaled), type = "scatter", name = 'v1/v4')

fig2 <- plot_ly(data = data.frame(reduced_M_scaled), x = ~v2, y = ~v4,
                text = rownames(reduced_M_scaled), type = "scatter", name = 'v2/v4')

fig3 <- plot_ly(data = data.frame(reduced_M_scaled), x = ~v3, y = ~v4,
                text = rownames(reduced_M_scaled), type = "scatter", name = 'v3/v4')



fig4 <- subplot(fig1, fig2, nrows = 1) 
p5 <- subplot(fig4, fig3, nrows = 2) 
p5
saveWidget(p5, "output/aygalic/PCA_axis_4.html", selfcontained = F, libdir = "lib")












##########################################################
############## PCA on the transposed matrix ##############
##########################################################
##########################################################

# We want to do PCA projecting the cell lines on a reduced gene space
tM_scaled = t(M_scaled)

cov2 = cov(tM_scaled)
decomp2 = eigen(cov2)
barplot(decomp2$values[1:50]/sum(decomp2$values), ylim = c(0,1))

# this is the projection space for the selected cancer type
P2 = as.matrix(decomp2$vectors[,1:2])






# next thing to do :
# Do PCA while hiding some of the classes from the data :
# only plot breast & lung or any combinaison

# plot all of them
project_cell_lines_for_cancer_type(space = P2)

#plot only breast cancer
project_cell_lines_for_cancer_type(c(7),space = P2)

#plot  breast cancer & a few other
project_cell_lines_for_cancer_type(c(11,19,15,7),space = P2)
project_cell_lines_for_cancer_type(c(13,15,7),space = P2)

project_cell_lines_for_cancer_type(c(8,9,7),space = P2)
project_cell_lines_for_cancer_type(c(11,10,7),space = P2)







