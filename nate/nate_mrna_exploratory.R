#REPLACE W YOUR PATH 
mrna = read.delim(file.path("../Dataset", "data_mrna_seq_rpkm.txt"), header = TRUE, comment.char = '#', nrows=100)
head(mrna)

#an idea for the data import/handling: c++ MPI parallel stuff :)) since its like ~53000x1100 samples 

#exploratory, re-formats column labels for cells into something easier for analyzing 
str_util <- function(str_arr){
  to_return = c(NULL)
  for(n in 1:length(str_arr)){
    temp <- unlist(strsplit(str_arr[n], "_"))
    trimmed = temp[2]
    if(length(temp)>=3){
      for(i in 3:length(temp)){
        trimmed<-paste(trimmed, temp[i],sep="_")
      }
    }
    to_return<-c(to_return,trimmed)
  }
  
  return(to_return);
}

#omit the hugo_symbol since it's not a cell 
trimmed = str_util(colnames(mrna)[-1])
#create a set 
cell_labels = trimmed[!duplicated(trimmed)]

#now let's make a factor and then check out the label counts 
freqs <- factor(trimmed, levels=cell_labels)
table(freqs)

#let's reduce the original dataframe into exclusively breast cells 
breast_labels = which(trimmed=="BREAST")+1 #+1 for the hugo_symbol...
breast_df = mrna[,breast_labels]

#rename row names since they were lost for some reason
rownames(breast_df)<-mrna[,1]

