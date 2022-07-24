library(dplyr)
library(ggplot2)

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
plot(log(row_var), pch = 16, cex = 0.3)
abline(h=log(threshold), col = "red")

data = data.frame(y = log(row_var), x = 1:length(row_var))


p <- ggplot(data = data, aes(x = x, y = y) ) + geom_point(size = .4) 
p <- p + geom_hline( yintercept = log(threshold), color = "red", size = 1)
#p <- p + ggtitle("Variability threshold")
p
p <- p + theme(
  panel.background = element_rect(fill = "transparent",
                                  colour = NA_character_), # necessary to avoid drawing panel outline
  panel.grid.major = element_blank(), # get rid of major grid
  panel.grid.minor = element_blank(), # get rid of minor grid
  plot.background = element_rect(fill = "transparent",
                                 colour = NA_character_), # necessary to avoid drawing plot outline
  legend.background = element_rect(fill = "transparent"),
  legend.box.background = element_rect(fill = "transparent"),
  legend.key = element_rect(fill = "transparent"),
  axis.text.x=element_blank(), #remove x axis labels
  axis.ticks.x=element_blank(), #remove x axis ticks
  axis.text.y=element_blank(),  #remove y axis labels
  axis.ticks.y=element_blank(),
  axis.title.x =element_blank(),
  axis.title.y =element_blank()
)

p
ggsave(
  plot = p,
  filename = "output/presentation/rpkm_threshold.png",
  bg = "transparent"
)


sum(row_var > threshold)


data_exp_var <- data_expression_clean[row_var > threshold,]


write.table(data_exp_var, file.path("Dataset", "1_rpkm.txt"),
            quote = FALSE, append = FALSE, sep = "\t", dec = ".", 
            row.names = TRUE, col.names = TRUE)


