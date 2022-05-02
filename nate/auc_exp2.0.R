source("./utils/nate_utils.R")
source("./utils/exp_base_script.R")

drugs.labelled<-read.csv("./some_stuff/drugs_labelled.csv")

#motivation: DIMENSIONALITY REDUCTION

#--------------------------------------
#VALIDATION OF PROPOSED STRAT - NOT NECESSARILY WHAT WE NEED TO JUSTIFY REDUCING TO MEAN 
#first check to see if data for a drug is guassian 
sort(table(cancer_freqs))

#lets check for LUNG and OVARY (corresponds to most samples, and cancer w mid range num samples)
lung.data<-block_dat("LUNG", auc)
lung.shapiro<-c(rep(NA), dim(lung.data)[1])
for(i in 1:dim(lung.data)[1]){
  sample<--as.numeric(lung.data[i,!is.na(lung.data[i,])])
  lung.shapiro[i]<- shapiro.test(sample)$p
}

x11()
par(mfrow=c(1,3))
hist(lung.shapiro, main="Histogram of p-values", xlab="auc")

max.index<-which(lung.shapiro == max(lung.shapiro))
min.index<-which(lung.shapiro == min(lung.shapiro))

hist(as.numeric(lung.data[max.index, !is.na(lung.data[max.index, ])]), main = rownames(lung.data)[max.index], xlab="auc")
hist(as.numeric(lung.data[min.index, !is.na(lung.data[min.index, ])]), main = rownames(lung.data)[min.index], xlab="auc")

dev.off()

#now repeat w ovarian cancer cells 

ovary.data<-block_dat("OVARY", auc)
ovary.shapiro<-c(rep(NA), dim(ovary.data)[1])
for(i in 1:dim(ovary.data)[1]){
  sample<--as.numeric(ovary.data[i,!is.na(ovary.data[i,])])
  ovary.shapiro[i]<- shapiro.test(sample)$p
}

x11()
par(mfrow=c(1,3))
hist(ovary.shapiro, main="Histogram of p-values")

max.index<-which(ovary.shapiro == max(ovary.shapiro))
min.index<-which(ovary.shapiro == min(ovary.shapiro))

hist(as.numeric(ovary.data[max.index, !is.na(ovary.data[max.index, ])]), main = rownames(ovary.data)[max.index])
hist(as.numeric(ovary.data[min.index, !is.na(ovary.data[min.index, ])]), main = rownames(ovary.data)[min.index])

dev.off()

#Note: in Ovarian and Pancreatic cases at least, most `gaussian` drug is of group: Approved 

#--------------------------------------
#Check to see if any relation with drug's approval status 

cancer<-"CENTRAL_NERVOUS_SYSTEM"

cancer.data<-block_dat(cancer, auc)
cancer.shapiro<-c(rep(NA), dim(cancer.data)[1])
for(i in 1:dim(cancer.data)[1]){
  sample<-as.numeric(cancer.data[i,!is.na(cancer.data[i,])])
  cancer.shapiro[i]<- shapiro.test(sample)$p
}

#min or max
shapiro.ordered<-sort(cancer.shapiro, decreasing=TRUE)

x11()
par(mfrow=c(4,4))

for(i in seq(1,16)){
  #20 bins
  breaks<-seq(0,1,1/100)
  index<-which(cancer.shapiro == shapiro.ordered[i])
  title<-paste(rownames(cancer.data)[index]," | pval=", format(cancer.shapiro[index], digits=2), " | ", drugs.labelled$Groups[index])
  hist(as.numeric(cancer.data[index, !is.na(cancer.data[index, ])]), main = title, xlab = "auc", breaks=breaks)
}

dev.off()


#--------------------------------------
#We might be more into reducing variance
#YEA THIS IS MORE WHAT WE WANT 

cancer<-"BREAST"

cancer.data<-block_dat(cancer, auc)
cancer.var<-c(rep(NA), dim(cancer.data)[1])
for(i in 1:dim(cancer.data)[1]){
  sample<-as.numeric(cancer.data[i,!is.na(cancer.data[i,])])
  cancer.var[i]<- var(sample)
}

#min or max
var.ordered<-sort(cancer.var, decreasing=TRUE)

x11()
par(mfrow=c(4,4))

for(i in seq(1,16)){
  #20 bins
  breaks<-seq(0,1,1/100)
  index<-which(cancer.var == var.ordered[i])
  title<-paste(rownames(cancer.data)[index]," | var=", format(cancer.var[index], digits=2), " | ", drugs.labelled$Groups[index])
  hist(as.numeric(cancer.data[index, !is.na(cancer.data[index, ])]), main = title, xlab = "auc", breaks=breaks)
}

dev.off()


#--------------------------------------

#we average a drug's performances over cancer types 

first<-block_dat(cancer_types[1], auc)
#won't affect the mean
for(i in 1:dim(first)[1]){
  row_dat <- first[i,]
  na_fill = NA
  if(length(row_dat[!is.na(row_dat)])>0){
    na_fill = mean(row_dat[!is.na(row_dat)])
  }
  first[i, is.na(first[i,])] = na_fill
}

training<-data.frame(rowMeans(first))
colnames(training)<-cancer_types[1];

for(c in cancer_types[-1]){
  print(c)
  cell_data<-block_dat(c, auc)
  
  #na fill
  for(i in 1:dim(cell_data)[1]){
    row_dat<-cell_data[i,]
    na_fill = NA
    if(length(row_dat[!is.na(row_dat)])>0){
      na_fill = mean(row_dat[!is.na(row_dat)])
    }
    cell_data[i, is.na(cell_data[i,])] = na_fill
  }
  
  averaged_col<-data.frame(rowMeans(cell_data));
  colnames(averaged_col)<-c;
  training<-cbind(training,averaged_col);
}

colnames(training)
unique(str_util(colnames(auc)))

View(block_dat("PROSTATE",auc))
