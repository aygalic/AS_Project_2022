breast_rpkm_data = read.delim(file.path("Dataset", "2_rpkm.txt"), header = TRUE, comment.char = '#')
breast_rpkm_data = t(breast_rpkm_data)
#colnames(breast_rpkm_data) = breast_rpkm_data[1,]
#breast_rpkm_data = breast_rpkm_data[-1,]

library(ElemStatLearn)
library(caret)
library(caTools)
library(e1071)

source("src/utilities.R")
source("robi/analysis_robi.rmd")

breast_auc_data = cancer_data_treatment_auc

#set the mean as value for NAs
rows = rowMeans(breast_auc_data, na.rm = TRUE)

for(i in 1:142)
{
  for(j in 1:43)
  {
    if(is.na(breast_auc_data[i,j])==TRUE)
    {
      breast_auc_data[i,j] = rows[i]
    }
  }
}
sum(is.na(breast_auc_data))
#no more NAs in the dataframe

#transpose to get cell lines in the rows
breast_auc_data = t(breast_auc_data)
#order alphabetically the cell lines
breast_auc_data = breast_auc_data[order(rownames(breast_auc_data)),]

M1 = as.matrix(breast_auc_data)
M2 = as.matrix(breast_rpkm_data)

# select row from M1 that are present in M2
indexes1 <- rownames(M1) %in% rownames(M2)
M1_ <- M1[indexes1,]
dim(M1_)

# select row from M2 that are present in M1
indexes2 <- rownames(M2) %in% rownames(M1)
M2_ <- M2[indexes2,]
dim(M2_)

result.k <- kmeans(M1_, centers=3)

reduced_M1_scaled = create_reduced_mat(M1_,2)


x11()
par(mfrow=c(1,2))
plot(M1_, col = result.k$cluster+1) #only the first 2 treatments
plot(reduced_M1_scaled, col=result.k$cluster+1, pch = 16, xlab = "PC1", ylab = "PC2") #on PC1/PC2 space


true_labels <- as.numeric(result.k$cluster)
M1_ = cbind(M1_, true_labels)

M2_ = cbind(M2_, true_labels)

new_1 = as.data.frame(breast_auc_data[indexes1,])
new_2 = as.data.frame(breast_rpkm_data[indexes2,])


# I TRY TO DO RANDOM FOREST but...
#inputData <- M2_[1:38, ] # training data
#testData <- M2_[15:42, ] # test data
#randomForest(formula =  true_labels ~ . , data = inputData, importance = TRUE, proximity = TRUE, ntree = 100, mtry = 6, plot = FALSE)
#-->too big data

result.k.rpkm <- kmeans(M2_, centers=3)
table(true = true_labels, pred = as.numeric(result.k.rpkm$cluster))
# clustering error between the 2 datasets
APER = (3+3+7+2)/(15+3+3+2+7+2+10)
APER


#remove columns (genes) with null values
ind = 0
for(i in 1:56318)
{
  if(new_2[,i] == rep(0,42))
    ind = cbind(ind,i)
}


new_2 = new_2[,-ind]
# from 56318 to 28953 genes 

# feature selection with thresholds (like LUCA)
# M <- scale(new_2) 

threshold_l <- 0.7 
threshold_h <- 100 


col_var = apply(new_2, 2, var) #apply over columns: variability along the genes (and clusters but that is into the bound interval)

plot((col_var))
abline(h=threshold_l, col='red')
abline(h=threshold_h, col='red')

sum(col_var > threshold_l & col_var < threshold_h)
plot(col_var[col_var > threshold_l & col_var < threshold_h])

new_2 <- new_2[ ,col_var > threshold_l & col_var < threshold_h]
# from 23953 to 12310 genes 


reduced_M2_scaled = create_reduced_mat(as.matrix(new_2),2)

x11()
par(mfrow=c(1,2))
plot(as.matrix(new_2), col = result.k$cluster+1, pch = 16) #only the first 2 treatments
plot(reduced_M2_scaled, col=result.k$cluster+1, pch = 16, xlab = "PC1", ylab = "PC2") #on PC1/PC2 space


# Classification
reduced_M2_scaled

dataset_reduced = data.frame(reduced_M2_scaled, cluster = factor(true_labels, levels = c(1,2,3)))



set.seed(123)
split = sample.split(dataset_reduced$cluster, SplitRatio = 0.8)

training_set = subset(dataset_reduced, split == TRUE)
test_set = subset(dataset_reduced, split == FALSE)

# Feature Scaling
training_set[,-3] = scale(training_set[,-3])
test_set[,-3] = scale(test_set[,-3])

library(e1071)
classifier = svm(formula = cluster ~ .,
                 data = training_set,
                 type = 'C-classification',
                 kernel = 'linear')

# Predicting the Test set results
y_pred = predict(classifier, newdata = test_set[-3])

# Making the Confusion Matrix
cm = table(test_set[,3], y_pred)
cm
#always classify as 1 (I also tried the non-linear case with kernel = 'radial' but the error is even higher!)

#misclassification error
(1+3)/8

# Plotting the training data set results

#X1  <- seq(min(training_set[,1]), max(training_set[,1]), length=200)
#X2  <- seq(min(training_set[,2]), max(training_set[,2]), length=200)

#grid_set <- expand.grid(X1,X2)
#colnames(grid_set) = c('v1', 'v2')
#y_grid = predict(classifier, newdata = grid_set)


#plot(training_set[,-3],
     #main = 'SVM (Training set)',
     #xlab = 'PC1', ylab = 'PC2',
     #xlim = range(X1), ylim = range(X2))

#contour(X1, X2, matrix(as.numeric(y_grid), length(X1), length(X2)), add = TRUE)

#points(grid_set, pch = '.', col = ifelse(y_grid == 1, 'coral1', 'aquamarine'))

#points(set, pch = 21, bg = ifelse(set[, 3] == 1, 'green4', 'red3'))



inTrain <- createDataPartition(y, p = .80, list = FALSE)[,1]

x_train <- x[ inTrain, ]
x_test  <- x[-inTrain, ]

y_train <- y[ inTrain]
y_test  <- y[-inTrain]


# Build a classifier from the auc clustering and train it on rpkm

## Support Vector Machine
dat <- data.frame(x = x_train, y=as.factor(y_train))
out <- svm(y~., data=dat , kernel ="linear",cost =10)
summary (out)

table(out$fitted , dat$y)

dat.te <- data.frame(x=x_test , y=as.factor(y_test))
pred.te <- predict (out , newdata =dat.te)
table(pred.te , dat.te$y)

## Random Forest

rpkm.rfe <- randomForest(cluster ~ ., 
                        data = dataset, 
                        importance = TRUE,
                        proximity = TRUE)

print(rpkm.rfe)
plot(rpkm)

