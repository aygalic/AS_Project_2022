#get our utils
source("./utils/nate_utils.R")

#pre-processing step/data import 
source("./utils/exp_base_script.R")

#----------------------------------------------------------------------------------------------
#from the imports 
cancer_types
sort(table(cancer_freqs), decreasing = TRUE)

querried_names = c("BREAST", "LARGE_INTESTINE", "CENTRAL_NERVOUS_SYSTEM", "SKIN", "OVARY", "PANCREAS")
#15 pairs 
name_pairs<-combn(querried_names,2)


#MOTIVATION: loop through pairs and plot drug pca on one chart
#----------------------------------------------------------------------------------------------
thresh = 0.04
which_cell<-function(data, thresh=0.05){
  max_val1 = max(data);
  data_copy = data[-(match(max_val1, data))]
  max_val2 = max(data_copy)
  
  if((max_val1-max_val2)<thresh){
    return("INDETERMINATE")
  }
  else{
    return (colnames(data)[match(max_val1,data)])
  }
}

#match the number of pairings for visualization
par(mfrow=c(5,3), mar=c(1,1,1,1))
for(i in 1:dim(name_pairs)[2]){
  #get pair 
  pair<-c(name_pairs[1,i], name_pairs[2,i])

  sub_auc<- block_dat(pair, auc)
  data_dim<-dim(sub_auc)
  
  #ALTERNATE 1 -> fill with const
  #na_fill = 0.5
  #sub_auc[is.na(sub_auc)]= na_fill
  
  #ALTERNATE 2 -> fill with column median
  for(i in 1:data_dim[2]){
    na_fill = median(sub_auc[!is.na(sub_auc[,i]), i])
    sub_auc[is.na(sub_auc[,i]),i] = na_fill
  }
  
  pc <- princomp(sub_auc,scores=T)
  
  #-------------------------------------------------------------
  #STUPID(?), CAUSE EVERY ITER WE CONTRADICT OURSELVES
  #this is only to generate plot colors 
  
  #classify what 'kind' of drug we have on a cancer-type basis 
  drugs<-data.frame(rowMeans(sub_auc[,grepl(pair[1], colnames(sub_auc))]))
  colnames(drugs)<-pair[1];
  
  for(c in pair[-1]){
    averaged_col = data.frame(rowMeans(sub_auc[,grepl(c, colnames(sub_auc))]));
    colnames(averaged_col)<-c;
    drugs<-cbind(drugs,averaged_col);
  }
  
  drug_max = c(NULL)
  for(i in 1:data_dim[1]){
    drug_max<-cbind(drug_max,which_cell(drugs[i,], thresh=thresh))
  }
  drug_max<-t(data.frame(drug_max));
  colnames(drug_max)<-"DRUG_MAX"
  drugs<-cbind(drugs,drug_max)
  
  #get some colors in here and assign them
  classes<-append(pair, "INDETERMINATE")
  col.cells <- rainbow(length(classes))
  col.labels <- rep(NA, data_dim[1])
  
  for(i in 1:length(col.labels)){
    col.labels[i] = col.cells[match(drugs[i,"DRUG_MAX"],classes)]
  }
  #-------------------------------------------------------------
  
  #just for the plot 
  if(pair[1]=="CENTRAL_NERVOUS_SYSTEM"){
    title<-paste("CNS", "-", pair[2])
  }
  else if(pair[2] == "CENTRAL_NERVOUS_SYSTEM"){
    title<-paste(pair[1], "-", "CNS")
  }
  else{
    title<-paste(pair[1], "-", pair[2])
  }
  
  #omit colors right now; col=col.labels
  plot(pc$scores[,1], pc$scores[,2],pch=16, col=col.labels, main=title, yaxt='n', xaxt='n', ylab = "", xlab="")
  
  #go 615x650 on dimension for plot save in R studio
}



