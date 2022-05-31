#PART 1: Basic plot of the features 

features = read.delim(file.path(getwd(), "physiological.csv"), header = TRUE, sep=",")
dim(features)
features[1,]

#NOTE: auc is average over not null drug-tested values

#replace empty string with NA 
features[features == ""] = NA

#encode the sexes 
features$SEX_ENCODED = as.numeric(factor(features$SEX)) - 1

#encode ethnicity 
features$ETHNICITY_ENCODED = as.numeric(factor(features$ETHNICITY)) - 1

#construct a matrix of plottables 
numerics = c("AGE", "DOUBLING_TIME", "MUTATION_RATE", "AUC", "SEX_ENCODED", "ETHNICITY_ENCODED")
features.numeric = features[, numerics]

x11()
plot(na.omit(features.numeric))



#------------------------------------------------------------
#PART 2: apply some color maps to scatter plots of interest
names(features)

#large data loss but still lots of samples 
interest = na.omit(features[,c("PATIENT_ID", "CANCER_TYPE","AGE", "AUC", "SEX", "ETHNICITY")])

dim(interest)[1]/dim(features)[1]

#generate some color maps based on ethnicity and sex combinations
table(interest$SEX)
table(interest$ETHNICITY)

Male = which(interest$SEX == "Male")
Female = which(interest$SEX == "Female")

AA = which(interest$ETHNICITY == "African_american")
A = which(interest$ETHNICITY == "Asian")
CA = which(interest$ETHNICITY == "Caucasian")

#colour map for sex
col1 = rep(NA,dim(interest)[1])
col1[Male] = 'red'
col1[Female] = 'blue'

#colour map for ethnicity
col2 = rep(NA,dim(interest)[1])
col2[AA] = 'black'
col2[A] = 'magenta'
col2[CA] = 'orange'

#plot age versus auc score 
x11()
par(mfrow = c(1,2))
plot(features$AGE, features$AUC, col = col1, pch = 19, main="Sex", xlab="Age", ylab="Average AUC")
plot(features$AGE, features$AUC, col = col2, pch = 19, main = "Ethnicity", xlab="Age", ylab="Average AUC")

#plot with all of them 
MaleAA = which(interest$SEX == "Male" & interest$ETHNICITY == "African_american")
FemaleAA = which(interest$SEX == "Female" & interest$ETHNICITY == "African_american")
MaleA = which(interest$SEX == "Male" & interest$ETHNICITY == "Asian")
FemaleA = which(interest$SEX == "Female" & interest$ETHNICITY == "Asian")
MaleCA = which(interest$SEX == "Male" & interest$ETHNICITY == "Caucasian")
FemaleCA = which(interest$SEX == "Female" & interest$ETHNICITY == "Caucasian")

col3 = rep(NA,dim(interest)[1])
col3[MaleAA] = 'green'
col3[FemaleAA] = 'red'
col3[MaleA] = 'blue'
col3[FemaleA] = 'orange'
col3[MaleCA] = 'pink'
col3[FemaleCA] = 'black'

shapes = rep(15,dim(interest)[1])
shapes[FemaleAA] = 16
shapes[MaleA] = 17
shapes[FemaleA] = 18
shapes[MaleCA] = 19
shapes[FemaleCA] = 20

x11()
plot(features$AGE, features$AUC, col = col3, pch = shapes, lwd=4, main="Sex and Ethnicity", xlab="Age", ylab="Average AUC")



#------------------------------------------------------------
#PART 3: restrict ourselves to a single cancer type 
interest.new = interest
interest.new$col1 = col1
interest.new$col2 = col2
interest.new$col3 = col3
interest.new$shapes = shapes

#build column with the string
source("../utils/nate_utils.R")
interest.new$CANCER = str_util(interest.new$PATIENT_ID)

#it turns out my method for parsing patient_id's yields same number of cancers as `CANCER_TYPE` field 
#so let's just use that (should be more accurate in theory)

sort(table(interest.new$CANCER_TYPE))

#pick a cancer 
cancer = c("Pancreatic Cancer")
cancer.data = interest.new[interest.new$CANCER_TYPE%in%cancer,]

table(cancer.data$SEX) 

#now lets plot 
x11()
par(mfrow=c(1,2))
plot(cancer.data$AGE, cancer.data$AUC, col = cancer.data$col1, main="Sex", xlab="Age", ylab="Average AUC", pch=19)
plot(cancer.data$AGE, cancer.data$AUC, col = cancer.data$col2, main="Ethnicity", xlab="Age", ylab="Average AUC", pch=19)


#------------------------------------------------------------
#PART 3.5: plot for top 16 cancers w same process
cancer.16 = names(sort(table(interest.new$CANCER_TYPE), decreasing = TRUE)[1:16])

x11()
par(mfrow(matrix(4,4)))
for(name in cancer.16){
  cancer.data = interest.new[interest.new$CANCER_TYPE%in%cancer,]
  plot(cancer.data$AGE, cancer.data$AUC, col = cancer.data$col1, main="Sex", xlab="Age", ylab="Average AUC", pch=19)
}
  



#------------------------------------------------------------
#PART4: color maps to the original pairwise plots 
#all features except cancer type, location, and id -> gonna use sex and ethnicity for colour maps 
interest2 = na.omit(features[,c("DOUBLING_TIME", "MUTATION_RATE","AGE", "AUC", "SEX", "ETHNICITY")])

#first we plot the with the sex map
sex_col = rep(NA,dim(interest2)[1])
Male = which(interest2$SEX == "Male")
Female = which(interest2$SEX == "Female")
sex_col[Male] = 'red'
sex_col[Female] = 'blue'

x11()
pairs(interest2[,c("DOUBLING_TIME","MUTATION_RATE","AGE","AUC")], col=sex_col, pch=19, main="Sex Mask")


#now let's do ethnicity
eth_col = rep('black',dim(interest2)[1])
AA = which(interest2$ETHNICITY == "African_american")
A = which(interest2$ETHNICITY == "Asian")
CA = which(interest2$ETHNICITY == "Caucasian")

eth_col[A] = 'red'
eth_col[CA] = 'pink'

x11()
pairs(interest2[,c("DOUBLING_TIME","MUTATION_RATE","AGE","AUC")], col=eth_col, pch=19, main="Ethnicity Mask")




#PART 4.5 log plots
#------------------------------------------------------------
#log plots
interest3 = interest2 
interest3$MUTATION_RATE = log(interest3$MUTATION_RATE)
interest3$DOUBLING_TIME = log(interest3$DOUBLING_TIME)


x11()
pairs(interest3[,c("DOUBLING_TIME","MUTATION_RATE","AGE","AUC")], col=sex_col, pch=19, main="Sex Mask, log")


x11()
pairs(interest3[,c("DOUBLING_TIME","MUTATION_RATE","AGE","AUC")], col=eth_col, pch=19, main="Ethnicity Mask, log", lwd=0.5)




#PART 5 Binning the ages and plotting w map based on that 
#------------------------------------------------------------



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

