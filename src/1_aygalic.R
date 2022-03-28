library(ggplot2)

library(plotly)  # fancy plots i guess

library(tidyr) # to use matrix for heatmap

setwd("~/OneDrive/polimi/COURSES/S8/APPLIED_STATS/AS_Project_2022")


##############################################
######### Code stolen from Roberta : #########
##############################################

#import data
data_patient= read.delim(file.path("Dataset", "data_clinical_patient.txt"), header = TRUE, comment.char = '#')

data_treatment_auc= read.delim(file.path("Dataset", 'data_drug_treatment_auc.txt'), header = TRUE, comment.char = '#')

data_sample= read.delim(file.path("Dataset",'data_clinical_sample.txt'), header = TRUE, comment.char = '#')

#select cells from breast carcinoma patients
selected_values = data_sample[data_sample$"CANCER_TYPE_DETAILED"=="Invasive Breast Carcinoma",]
selected_cells = selected_values$SAMPLE_ID

selected_cells = na.omit(selected_cells) # I get 57 values

indexes = match(selected_cells, colnames(data_treatment_auc))

sum(is.na(indexes)) 
#there are 12 NA so we end up with data on AUC for 45 cell lines with breast cancer

cancer_data_treatment_auc = data_treatment_auc[, c(1, na.omit(indexes)) ]

perc_NA_row = rowMeans(is.na(cancer_data_treatment_auc))*100
#10% tolerance on the rows and I see if the mean percentage is 
max_perc_row = max(perc_NA_row)
while(max_perc_row > 10 && mean(perc_NA_row) > 2)
{
  index = match(max_perc_row, perc_NA_row)
  cancer_data_treatment_auc = cancer_data_treatment_auc[-c(index),]
  perc_NA_row = perc_NA_row[-c(index)]
  max_perc_row = max(perc_NA_row)
}


#update the percentages for columns
perc_NA_col = colMeans(is.na(cancer_data_treatment_auc))*100


###############################################
############### Original Code : ###############
###############################################


#standardization ???
auc_data_standardized <- scale(cancer_data_treatment_auc[-1,-1])

# a few plot : std vs non std
boxplot(auc_data_standardized, use.cols = TRUE)
boxplot(cancer_data_treatment_auc[-1,-1], use.cols = TRUE)




data_matrix <- as.matrix(cancer_data_treatment_auc[-1,-1])
#replace NAs with .5
data_matrix[is.na(data_matrix)] <- 0.5

fig <- plot_ly(z = data_matrix, type = "heatmap")
fig


# with normalized data
data_matrix_std <- scale(as.matrix(cancer_data_treatment_auc[-1,-1]))
data_matrix_std[is.na(data_matrix_std)] <- 0.5

fig_std <- plot_ly(z = data_matrix_std, type = "heatmap")
fig_std
