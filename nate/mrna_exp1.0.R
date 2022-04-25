#REPLACE W YOUR PATH 
mrna = read.delim(file.path("../Dataset", "data_mrna_seq_rpkm.txt"), header = TRUE, comment.char = '#', nrows=1000)


#an idea for the data import/handling: c++ MPI parallel stuff :)) since its like ~53000x1100 samples 

#exploratory, re-formats column labels for cells into something easier for analyzing 
str_util <- function(str_arr){
  to_return = c(NULL)
  for(n in 1:length(str_arr)){
    temp <- unlist(strsplit(str_arr[n], "_"))
    if(length(temp)>1){
      trimmed = temp[2]
      if(length(temp)>=3){
        for(i in 3:length(temp)){
          trimmed<-paste(trimmed, temp[i],sep="_")
        }
      }
      to_return<-c(to_return,trimmed)
    }
    else{
      to_return = temp[1]
    }
  }
  return(to_return)
}

#----------------------------------------------------------------------------------------------
#EXPLORATORY 

#omit the hugo_symbol since it's not a cell 
trimmed = str_util(colnames(mrna)[-1])
#create a set 
cell_labels = unique(trimmed)

#now let's make a factor and then check out the label counts 
freqs <- factor(trimmed, levels=cell_labels)
freq_table<-table(freqs)
red_table<-freq_table[which(freq_table>20)]

pie(red_table, col = hcl.colors(length(red_table), "BluYl"), 
    main = "Breakdown of cancers present in dataset")

#let's reduce the original dataframe into exclusively breast cells 
breast_labels = which(trimmed=="BREAST")+1 #+1 for the hugo_symbol...
breast_df = mrna[,breast_labels]

#rename row names since they were lost for some reason
rownames(breast_df)<-mrna[,1]

#----------------------------------------------------------------------------------------------
#LET'S TRY TO EXPEDIATE THE PROCESS W THIS FUNCTION FROM OTHER FILE 
block_dat <-function(cancer_types_, cancer_df){
  #first pass 
  type = cancer_types_[1]
  selected_cells = colnames(cancer_df)[grepl(type, colnames(cancer_df))]
  cell_indices = match(selected_cells,colnames(cancer_df))
  to_return = data.frame(cancer_df[, cell_indices])
  if(length(selected_cells)==1){
    colnames(sub_auc)<-selected_cells
  }
  
  for(i in 2:length(cancer_types_)){
    type = cancer_types_[i]
    selected_cells = colnames(cancer_df)[grepl(type, colnames(cancer_df))]
    cell_indices = match(selected_cells,colnames(cancer_df))
    sub_auc<-data.frame(cancer_df[, cell_indices])
    if(length(selected_cells)==1){
      colnames(sub_auc)<-selected_cells
    }
    to_return<-cbind(to_return, sub_auc)
  }
  
  return(to_return)
}

mrna_sub = block_dat(c("BONE", "BREAST"), mrna)

