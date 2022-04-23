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


block_dat <-function(cancer_types_, cancer_df){
  #first pass 
  type = cancer_types_[1]
  selected_cells = colnames(cancer_df)[grepl(type, colnames(cancer_df))]
  cell_indices = match(selected_cells,colnames(cancer_df))
  to_return = data.frame(cancer_df[, cell_indices])
  if(length(selected_cells)==1){
    colnames(sub_auc)<-selected_cells
  }
  
  if(length(cancer_types_)>1){
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
  }
  
  return(to_return)
}


  
