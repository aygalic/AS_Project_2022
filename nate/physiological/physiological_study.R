#PART 0: Check some properties of the data 
features = read.delim(file.path(getwd(), "physiological.csv"), header = TRUE, sep=",")
dim(features)
features[1,]

#check what the losses would be from dropping NA's for specific features
sum(is.na(features$DOUBLING_TIME))/length(features$DOUBLING_TIME) #40% loss

sum(is.na(features$MUTATION_RATE))/length(features$MUTATION_RATE) #65% loss

sum(is.na(features$ETHNICITY))/length(features$ETHNICITY) #0% loss

sum(is.na(features$AGE))/length(features$AGE) #55% loss


#PART 1: Basic plot of the features 
#replace empty string with NA 
features[features == ""] = NA

#encode the sexes 
features$SEX_ENCODED = as.numeric(factor(features$SEX)) - 1

#encode ethnicity 
features$ETHNICITY_ENCODED = as.numeric(factor(features$ETHNICITY)) - 1

#bin the ages
library(dplyr)
features$AGE_ENCODED = as.numeric(factor(cut(features$AGE, breaks = c(0,25,50,75,100))))-1


#construct a matrix of plottables 
numerics = c("AGE", "DOUBLING_TIME", "MUTATION_RATE", "AUC")
features.numeric = features[, numerics]

x11()
plot(na.omit(features.numeric))



#------------------------------------------------------------
#PART 2: apply some color maps to scatter plots of interest
names(features)

#large data loss but still lots of samples.  DEGREE TO WHICH LOSSES OCCUR QUANTIFIED IN PART 0
#got rid of patient id and the primary site for now 
interest = na.omit(features[,-c(1,8)])

dim(interest)[1]/dim(features)[1]

#generate some color maps based on ethnicity and sex combinations
table(interest$SEX)
table(interest$ETHNICITY)
table(interest$AGE_ENCODED)

Male = which(interest$SEX == "Male")
Female = which(interest$SEX == "Female")

AA = which(interest$ETHNICITY == "African_american")
A = which(interest$ETHNICITY == "Asian")
CA = which(interest$ETHNICITY == "Caucasian")

B0 = which(interest$AGE_ENCODED == 0)
B1 = which(interest$AGE_ENCODED == 1)
B2 = which(interest$AGE_ENCODED == 2)
B3 = which(interest$AGE_ENCODED == 3)

#colour map for sex
col_sex = rep(NA,dim(interest)[1])
col_sex[Male] = 'red'
col_sex[Female] = 'blue'

#colour map for ethnicity
col_eth = rep(NA,dim(interest)[1])
col_eth[AA] = 'black'
col_eth[A] = 'red'
col_eth[CA] = 'pink'

#colour map for age 
col_age = rep(NA,dim(interest)[1])
col_age[B0] = 'blue'
col_age[B1] = 'red'
col_age[B2] = 'pink'
col_age[B3] = 'black'

#plot some maps 
x11()
pairs(interest[,c("DOUBLING_TIME","MUTATION_RATE","AUC", "AGE")], col=col_sex, pch=19, main="Sex Mask")

#plot some maps 
x11()
pairs(interest[,c("DOUBLING_TIME","MUTATION_RATE","AUC", "AGE")], col=col_eth, pch=19, main="Ethnicity Mask")

#plot some maps 
x11()
pairs(interest[,c("DOUBLING_TIME","MUTATION_RATE","AUC")], col=col_age, pch=19, main="Age Mask")


#------------------------------------------------------------
#PART 3: restrict ourselves to a single cancer type

#save colours for easy access
interest$col_sex = col_sex
interest$col_eth = col_eth
interest$col_age = col_age

sort(table(interest$CANCER_TYPE))

#pick a cancer 
cancer = c("Pancreatic Cancer")
cancer.data = interest[interest$CANCER_TYPE%in%cancer,]

table(cancer.data$SEX) 

#now lets plot 
x11()
par(mfrow=c(1,2))
plot(cancer.data$AGE, cancer.data$AUC, col = cancer.data$col_sex, main="Sex", xlab="Age", ylab="Average AUC", pch=19, lwd=2)
plot(cancer.data$AGE, cancer.data$AUC, col = cancer.data$col_eth, main="Ethnicity", xlab="Age", ylab="Average AUC", pch=19, lwd=2)


#------------------------------------------------------------
#PART 3.5: plot for top n^2 cancers w same process
size = 9
cancers_chosen = names(sort(table(interest$CANCER_TYPE), decreasing = TRUE)[1:size])

#age dependence plots 
x11()
par(mfrow=c(sqrt(size),sqrt(size)))
for(name in cancers_chosen){
  cancer.data = interest[interest$CANCER_TYPE%in%name,]
  plot(cancer.data$AGE, cancer.data$AUC, col = cancer.data$col_sex, main=name, xlab="Age", ylab="Average AUC", pch=19, lwd=2)
}

x11()
par(mfrow=c(sqrt(size),sqrt(size)))
for(name in cancers_chosen){
  cancer.data = interest[interest$CANCER_TYPE%in%name,]
  plot(cancer.data$AGE, cancer.data$AUC, col = cancer.data$col_eth, main=name, xlab="Age", ylab="Average AUC", pch=19, lwd=2)
}


#NOTE: the choice to plot auc vs age w color maps suggests you think theres a relation between them
# this would mean you get a multiple regression model for auc vs age w groups designated by colour var
# if auc is the target we should be plotting (mut.rate, doub.time, age) x (auc) x {col.maps} number 
# of plots to `test` all our theories 

#wait till we log transform tho 


#PART 4 log plots
#------------------------------------------------------------
#log plots
interest.log = interest
interest.log$MUTATION_RATE = log(interest.log$MUTATION_RATE)
interest.log$DOUBLING_TIME = log(interest.log$DOUBLING_TIME)


x11()
pairs(interest.log[,c("DOUBLING_TIME","MUTATION_RATE","AGE","AUC")], col=interest$col_sex, pch=19, main="Sex Mask, log")

x11()
pairs(interest.log[,c("DOUBLING_TIME","MUTATION_RATE","AGE","AUC")], col=interest$col_eth, pch=19, main="Ethnicity Mask, log")

x11()
pairs(interest.log[,c("DOUBLING_TIME","MUTATION_RATE","AUC")], col=interest$col_age, pch=19, main="Ethnicity Mask, log")



#top cancers again 
#doubling time 
x11()
par(mfrow=c(sqrt(size),sqrt(size)))
for(name in cancers_chosen){
  cancer.data = interest.log[interest.log$CANCER_TYPE%in%name,]
  plot(cancer.data$DOUBLING_TIME, cancer.data$AUC, col = cancer.data$col_sex, main=name, xlab="log doubling time", ylab="Average AUC", pch=19, lwd=2)
}

x11()
par(mfrow=c(sqrt(size),sqrt(size)))
for(name in cancers_chosen){
  cancer.data = interest.log[interest.log$CANCER_TYPE%in%name,]
  plot(cancer.data$DOUBLING_TIME, cancer.data$AUC, col = cancer.data$col_eth, main=name, xlab="log doubling time", ylab="Average AUC", pch=19, lwd=2)
}

x11()
par(mfrow=c(sqrt(size),sqrt(size)))
for(name in cancers_chosen){
  cancer.data = interest.log[interest.log$CANCER_TYPE%in%name,]
  plot(cancer.data$DOUBLING_TIME, cancer.data$AUC, col = cancer.data$col_age, main=name, xlab="log doubling time", ylab="Average AUC", pch=19, lwd=2)
}



#mutation rate 
x11()
par(mfrow=c(sqrt(size),sqrt(size)))
for(name in cancers_chosen){
  cancer.data = interest.log[interest.log$CANCER_TYPE%in%name,]
  plot(cancer.data$MUTATION_RATE, cancer.data$AUC, col = cancer.data$col_sex, main=name, xlab="log mutation rate", ylab="Average AUC", pch=19, lwd=2)
}

x11()
par(mfrow=c(sqrt(size),sqrt(size)))
for(name in cancers_chosen){
  cancer.data = interest.log[interest.log$CANCER_TYPE%in%name,]
  plot(cancer.data$MUTATION_RATE, cancer.data$AUC, col = cancer.data$col_eth, main=name, xlab="log mutation rate", ylab="Average AUC", pch=19, lwd=2)
}

x11()
par(mfrow=c(sqrt(size),sqrt(size)))
for(name in cancers_chosen){
  cancer.data = interest.log[interest.log$CANCER_TYPE%in%name,]
  plot(cancer.data$MUTATION_RATE, cancer.data$AUC, col = cancer.data$col_age, main=name, xlab="log mutation rate", ylab="Average AUC", pch=19, lwd=2)
}





#PART 6: added specificity for the auc averaging 
############################################################################
######### SECTION FOR REDEFINING AUC AS AVERAGE OVER SPECIFIC DRUGS ########
#patient ids which are compatible for referencing in auc dataframe 
cancer.ids = cancer.data$PATIENT_ID

source("utils/exp_base_script.R")

#not all the cells in our patient dataset have been drug-tested, 
cancer.ids_new = cancer.ids[cancer.ids%in%colnames(auc)]
length(cancer.ids_new)/length(cancer.ids)

cancer.auc_data = auc[,cancer.ids_new]

#reduce original cancer.data to the new rows 
cancer.data_new = cancer.data[cancer.data$PATIENT_ID%in%cancer.ids_new,]

#################################################
# Aside: this `cancer.auc_data` is extremely similar to the result
# of calling my `block_dat` function with cancer_type = "LARGE_INTESTINE".
# here we've just modified the paradigm to use the dataset's original 
# labels for cancers instead of parsing the ids...
#################################################
#pf
cancer.ids%in%colnames(block_dat("LARGE_INTESTINE", auc))

#rewrite the AUC column to be the average over our new data 
#THIS IS THE PART WHERE YOU DEFINE THE NEW AUC AVERAGE E.G. REDUCE ROWS TO SPECIFIC
#DRUGS OR SOMETHING 
cancer.data_new$AUC = colMeans(cancer.auc_data, na.rm=TRUE)

############################################################################
############################################################################

