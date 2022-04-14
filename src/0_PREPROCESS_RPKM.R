library(dplyr)

setwd("~/OneDrive/polimi/COURSES/S8/APPLIED_STATS/AS_Project_2022")

data_expression= read.delim(file.path("Dataset", "data_mrna_seq_rpkm.txt"), header = TRUE, comment.char = '#')


# PROBLEM 1: DUPLICATES 

# Keep only non duplicated element
#
# We still end up with a few duplicate element especially Y-RNA
# But we decided not to deal with this for now
duplicated_values <- duplicated(data_expression)
data_expression_clean <- distinct(data_expression, .keep_all=TRUE)


write.table(data_expression_clean, file.path("Dataset", "0_rpkm.txt"),
            quote = FALSE, append = FALSE, sep = "\t", dec = ".", 
            row.names = TRUE, col.names = TRUE)

# PROBLEM 2: Removing lines with low variance

# we first create a separated matrix of normalized data
M <- scale(data_expression_clean[,-1])


threshold <- 0.5
row_var = apply(M, 1, var)
plot(row_var)

sum(row_var > threshold)
plot(row_var[row_var > threshold])

data_exp_var <- data_expression_clean[row_var > threshold,]


write.table(data_exp_var, file.path("Dataset", "1_rpkm.txt"),
            quote = FALSE, append = FALSE, sep = "\t", dec = ".", 
            row.names = TRUE, col.names = TRUE)


