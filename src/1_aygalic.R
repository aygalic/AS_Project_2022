library(ggplot2)
library(plotly)  # interactive plots 

setwd("~/OneDrive/polimi/COURSES/S8/APPLIED_STATS/AS_Project_2022")


##############################################
######### Code stolen from Roberta : #########
##############################################

#import data
data_patient= read.delim(file.path("Dataset", "data_clinical_patient.txt"), header = TRUE, comment.char = '#')

data_sample= read.delim(file.path("Dataset",'data_clinical_sample.txt'), header = TRUE, comment.char = '#')

original_data_mrna = read.delim(file.path("Dataset", "data_mrna_seq_rpkm.txt"), header = TRUE, comment.char = '#', nrows=100)


#select cells from breast carcinoma patients
selected_values = data_sample[data_sample$"CANCER_TYPE_DETAILED"=="Invasive Breast Carcinoma",]
selected_cells = selected_values$SAMPLE_ID
cat(length(selected_cells), "cells selected from cancer type Invasive Breast Carcinoma")
selected_cells = na.omit(selected_cells)
cat(length(selected_cells), "cells selected after omiting NAs")

indexes_mrna = match(selected_cells, colnames(original_data_mrna))
mrna_data = original_data_mrna[, c(1, na.omit(indexes_mrna)) ]


###############################################
############### Original Code : ###############
###############################################

cancer_data_mrna <- mrna_data[,-1]
#standardization ???
cancer_data_mrna_std <- scale(cancer_data_mrna)


# The data I was actually supposed to plot :
boxplot(cancer_data_mrna, use.cols = TRUE)
boxplot(cancer_data_mrna_std, use.cols = TRUE)

Y = mrna_data[,1]
X = colnames(cancer_data_mrna)

# with normalized data 
fig_std_mrna <- plot_ly(x = X, y = Y, z = cancer_data_mrna_std, type = "heatmap")
fig_std_mrna

# with non normalized data 
fig_mrna <- plot_ly(x = X, y = Y, z = as.matrix(cancer_data_mrna), type = "heatmap")
fig_mrna
