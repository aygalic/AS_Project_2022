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

data_patient = read.delim(file.path("../Dataset", "data_clinical_patient.txt"), header = TRUE, comment.char = '#')
data_sample = read.delim(file.path("../Dataset", "data_clinical_sample.txt"), header = TRUE, comment.char = '#')
data_treatment_auc = read.delim(file.path("../Dataset", 'data_drug_treatment_auc.txt'), header = TRUE, comment.char = '#')


#--------------------------------------------------------------------------------------------
#INTERSECTION OF DATASETS BASED ON PATIENT ID//PREPROCESSING  
ids <-na.omit(data_patient$PATIENT_ID)
ids

#note: we should intersect all datasets we plan to use, i.e. also consider data_mrna_seq_rpkm.txt

#auc 
overlap_indices_auc <- match(ids, colnames(data_treatment_auc))
overlap_indices_auc<-na.omit(overlap_indices_auc)
reduced_auc <-data_treatment_auc[,c(overlap_indices_auc)]
rownames(reduced_auc)<-data_treatment_auc$ENTITY_STABLE_ID

#samples (for some reason)
overlap_indices_sample <- match(ids, data_sample$PATIENT_ID)
overlap_indices_sample<-na.omit(overlap_indices_sample)
reduced_sample <-data_sample[overlap_indices_sample,]

#unique names from reduced 
reduced_ids = colnames(reduced_auc)
cancer_types = unique(str_util(reduced_ids))

#choose cancers of interest, maybe top 5-10 OPTIONAL 
freqs <- factor(str_util(reduced_ids), levels=cancer_types)
table(freqs)

#----------------------------------------------------------------------------------------------
#SOME CODE I WROTE WHICH IS IMPLEMENTED IN A LOOP BELOW 

type = "SALIVARY_GLAND"
#selected_cells = data_patient$PATIENT_ID[grepl(type,data_patient$PATIENT_ID)]
selected_cells = colnames(reduced_auc)[grepl(type, colnames(reduced_auc))]
#cell_indices = na.omit(match(selected_cells, colnames(reduced_auc)))
cell_indices = match(selected_cells,colnames(reduced_auc))
sub_auc= data.frame(reduced_auc[, cell_indices])
print(unique(str_util(colnames(sub_auc))))

#----------------------------------------------------------------------------------------------

#ORGANIZE THE MAIN DATASET INTO BLOCKS BASED ON CANCER TYPES 

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

querried_names = c("BREAST", "BONE", "LIVER")
sub_auc<- block_dat(querried_names, reduced_auc)



