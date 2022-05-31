#extract rpkm for breast cancer
source("nate/utils/nate_utils.R")
library(dplyr)


data_expression= read.delim(file.path("Dataset", "data_mrna_seq_rpkm.txt"), header = TRUE, comment.char = '#')


# PROBLEM 1: DUPLICATES 

# Keep only non duplicated element
#
# We still end up with a few duplicate element especially Y-RNA
# But we decided not to deal with this for now
data_rpkm = data_expression[,1]
breast_canc_cells = colnames(data_expression)[1]
for(i in 1:length(colnames(data_expression)))
{
  if(str_util(colnames(data_expression)[i]) == "BREAST")
  {
    data_rpkm = cbind(data_rpkm, data_expression[,i])
    breast_canc_cells = cbind(breast_canc_cells, colnames(data_expression)[i])
  }
}
colnames(data_rpkm) = breast_canc_cells
row.names(data_rpkm) <- data_rpkm[,1]
data_rpkm = data_rpkm[,-1]

#deal with duplicates
data_rpkm = as.data.frame(data_rpkm)
duplicated_values <- duplicated(data_rpkm)

data_expression_clean <- distinct(data_rpkm, .keep_all=TRUE)

write.table(data_expression_clean, file.path("Dataset", "0_rpkm.txt"),
            quote = FALSE, append = FALSE, sep = "\t", dec = ".", 
            row.names = TRUE, col.names = TRUE)

# PROBLEM 2: Removing lines with low variance

# we first create a separated matrix of normalized data
data_expression_clean = data.matrix(data_expression_clean)
M <- scale((data_expression_clean[,-1]))

threshold <- 0.5
row_var = apply(M, 1, var)
plot(row_var)

sum(row_var > threshold)
plot(row_var[row_var > threshold])

data_exp_var <- data_expression_clean[row_var > threshold,]

#select cancer cells



write.table(data_rpkm, file.path("Dataset", "2_rpkm.txt"),
            quote = FALSE, append = FALSE, sep = "\t", dec = ".", 
            row.names = TRUE, col.names = TRUE)


