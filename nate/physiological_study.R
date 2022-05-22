features = read.delim(file.path(getwd(), "physiological.csv"), header = TRUE, sep=",")
dim(features)
features[1,]

#encode the sexes 
features$SEX_ENCODED = as.numeric(factor(features$SEX)) - 1

#encode ethnicity 
features$ETHNICITY_ENCODED = as.numeric(factor(features$ETHNICITY)) - 1


#construct a matrix of plottables 
numerics = c("AGE", "DOUBLING_TIME", "MUTATION_RATE", "AUC", "SEX_ENCODED", "ETHNICITY_ENCODED")
features.numeric = features[, numerics]

x11()
plot(na.omit(features.numeric))

#no immediately apparent trends between variables i've selected, then again the auc is an average over all 
#drugs which might not be a great indication of performance.  could consider max auc and other stuff, or 
#a reduced subset of tested drugs 
