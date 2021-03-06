library(ggplot2)
library(plotly)  # interactive plots 

setwd("~/OneDrive/polimi/COURSES/S8/APPLIED_STATS/AS_Project_2022")
# all the useful functions for heatmaps
source("src/utilities.R")

#import data
data_patient_= read.delim(file.path("Dataset", "data_clinical_patient.txt"), header = TRUE, comment.char = '#')
data_sample_= read.delim(file.path("Dataset",'data_clinical_sample.txt'), header = TRUE, comment.char = '#')
original_data_mrna_ = read.delim(file.path("Dataset", "data_mrna_seq_rpkm.txt"), header = TRUE, comment.char = '#', nrows=5000)
original_data_mrna_ = read.delim(file.path("Dataset", "1_rpkm.txt"), header = TRUE, comment.char = '#', nrows=5000)


# Get a table of all cancer types and occurence
cancer_types <- as.data.frame(table(data_sample_$"CANCER_TYPE_DETAILED"), stringsAsFactors = FALSE)
cancer_types <- cancer_types[order(cancer_types$Freq, decreasing = TRUE),]
names(cancer_types)<-c("Factor", "Freq")



###########################################
############ Building Heatmaps ############ 
###########################################

# creating a matrix will all cancer types
#
# We could also create a matrix with specified cancer types only
cancer_type_selection <- c(1:(length(cancer_types$Factor)-1))
M = Build_matrix_for_multiple_cancer_types(cancer_type_selection)
tags = M$Tags
M <- M$Mat
M_scaled <- as.matrix(scale(M))


# plot heatmap
p1 <- plot_ly(x=colnames(M_scaled), y=rownames(M_scaled), z = M_scaled, type = "heatmap")
p1

saveWidget(p1, "output/aygalic/heatmap.html", selfcontained = F, libdir = "lib")


# Dendrogram
library(heatmaply)
heatmaply(M_scaled, k_col = length(cancer_types$Factor) -1)

# breaks the display a little, this is the fix :
dev.off()



