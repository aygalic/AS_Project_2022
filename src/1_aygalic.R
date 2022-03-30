library(ggplot2)
library(plotly)  # interactive plots 

setwd("~/OneDrive/polimi/COURSES/S8/APPLIED_STATS/AS_Project_2022")


#import data
data_patient_= read.delim(file.path("Dataset", "data_clinical_patient.txt"), header = TRUE, comment.char = '#')
data_sample_= read.delim(file.path("Dataset",'data_clinical_sample.txt'), header = TRUE, comment.char = '#')
original_data_mrna_ = read.delim(file.path("Dataset", "data_mrna_seq_rpkm.txt"), header = TRUE, comment.char = '#', nrows=100)


plot_heatmap <- function(cancer_name,
                         data_patient = data_patient_,
                         data_sample = data_sample_,
                         original_data_mrna = original_data_mrna_,
                         std=1){
  
  #select cells from breast carcinoma patients
  selected_values = data_sample[data_sample$"CANCER_TYPE_DETAILED"==cancer_name,]
  selected_cells = selected_values$SAMPLE_ID
  selected_cells = na.omit(selected_cells)
  indexes_mrna = match(selected_cells, colnames(original_data_mrna))
  mrna_data = original_data_mrna[, c(1, na.omit(indexes_mrna)) ]
  
  #select values only
  data <- as.matrix(mrna_data[,-1])
  
  if(std){data <- scale(data)}
  
  Y = mrna_data[,1]
  X = colnames(data)
  
  fig <- plot_ly(x = X, y = Y, z = data, type = "heatmap")%>%
    layout(title = cancer_name)
  fig
}




cancer_types <- as.data.frame(table(data_sample$"CANCER_TYPE_DETAILED"), stringsAsFactors = FALSE)
cancer_types <- cancer_types[order(cancer_types$Freq, decreasing = TRUE),]
names(cancer_types)<-c("Factor", "Freq")

plot_heatmap(cancer_name = "Invasive Breast Carcinoma")
plot_heatmap(cancer_name = "Small Cell Lung Cancer")

plot_heatmap(cancer_name = cancer_types$Factor[1])
plot_heatmap(cancer_name = cancer_types$Factor[2])
plot_heatmap(cancer_name = cancer_types$Factor[3])
plot_heatmap(cancer_name = cancer_types$Factor[4])
plot_heatmap(cancer_name = cancer_types$Factor[5])

