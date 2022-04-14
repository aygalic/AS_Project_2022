library(dplyr)

setwd("~/OneDrive/polimi/COURSES/S8/APPLIED_STATS/AS_Project_2022")

data_expression= read.delim(file.path("Dataset", "data_mrna_seq_rpkm.txt"), header = TRUE, comment.char = '#')


## PROBLEM 1: DUPLICATES 


duplicated_values <- duplicated(data_expression)
sum(duplicated_values)


data_expression[duplicated_values,]


data_expression_clean <- distinct(data_expression, .keep_all=TRUE)


write.table(data_expression_clean, file.path("Dataset", "0_rpkm.txt"),
            quote = FALSE, append = FALSE, sep = "\t", dec = ".", 
            row.names = TRUE, col.names = TRUE)

