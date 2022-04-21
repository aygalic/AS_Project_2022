#imports 
source("./utils/nate_utils.R")
source("./utils/exp_base_script.R")

#now import custom tki_indicator from tki.ipynb
path<-"/Users/nathanielnethercott/Desktop/School/Polimi/2021:2022/AS/AS_Project_2022/nate/some_stuff"
tki<-read.delim(file.path(path, "tki_labels.txt"), header = FALSE, comment.char = '#')
colnames(tki)<-"tki"

#same idea as auc_exp1.1 but we have a new coloring scheme!! 

querried_names = c("BREAST", "LARGE_INTESTINE", "CENTRAL_NERVOUS_SYSTEM", "SKIN", "OVARY", "PANCREAS")
querried_names = unique(str_util(colnames(auc)))

#15 pairs 
name_pairs<-combn(querried_names,2)


#match the number of pairings for visualization
par(mfrow=c(5,3), mar=c(1,1,1,1))
for(i in 1:dim(name_pairs)[2]){
  #get pair 
  pair<-c(name_pairs[1,i], name_pairs[2,i])
  
  pair<-querried_names
  
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
  
  #PCA 
  pc <- princomp(sub_auc,scores=T)
  
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
  
  col.labels<-rep('green',dim(tki)[1])
  for(i in 1:dim(tki)[1]){
    if(tki[i,] == 1){
      col.labels[i] = 'red'
    }
  }
  
  #omit colors right now; col=col.labels
  plot(pc$scores[,1], pc$scores[,2],pch=16, col=col.labels)
  
  #go 615x650 on dimension for plot save in R studio
}



