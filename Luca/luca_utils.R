#source(file.path("nate", "utils", "nate_utils.R"))
#import nate's function
#setwd("/Users/lucamainini/Documents/GitHub/AS_Project_2022")


import_dataset <- function(path){
  #' Import AUC and RPKM database intersecting them
  #' 
  #' @param path of the folder with dataframe
  #' @return a list with AUC and RPKM data frame
  #' @examples
  #' path = "Dataset"
 
  #DATA SET IMPORT
  data_patient = read.delim(file.path(path, "data_clinical_patient.txt"), header = TRUE, comment.char = '#')
  data_sample = read.delim(file.path(path, "data_clinical_sample.txt"), header = TRUE, comment.char = '#')
  data_treatment_auc = read.delim(file.path(path, "data_drug_treatment_auc.txt"), header = TRUE, comment.char = '#')
  data_rpkm = read.delim(file.path(path, "data_mrna_seq_rpkm.txt"), header = TRUE, comment.char = '#')
  
  #--------------------------------------------------------------------------------------------
  #INTERSECTION OF DATASETS BASED ON PATIENT ID/PREPROCESSING  
  #THIS IS ONE POSSIBLE PATH 
  
  #get names of patient-ids for cross-comparison
  ids <-na.omit(data_patient$PATIENT_ID)
  
  #we should intersect all datasets we plan to use on the basis of patient sample-id
  
  #intersection, carry-forward 
  auc_ids <- match(ids, colnames(data_treatment_auc))
  auc_ids<-na.omit(auc_ids)
  reduced_auc_ids <-colnames(data_treatment_auc)[auc_ids]
  
  sample_ids <- match(reduced_auc_ids, data_sample$PATIENT_ID)
  sample_ids<-na.omit(sample_ids)
  reduced_sample_ids <-data_sample$PATIENT_ID[sample_ids]
  
  mrna_ids<-match(reduced_sample_ids, colnames(data_rpkm))
  mrna_ids<-na.omit(mrna_ids)
  reduced_mrna_ids <-colnames(data_rpkm)[mrna_ids]
  
  reduced_ids<-reduced_mrna_ids
  
  #REASSIGNMENT 
  auc <-data_treatment_auc[,c(reduced_ids)]
  rownames(auc)<-data_treatment_auc$ENTITY_STABLE_ID
  
  rpkm <-data_rpkm[,c(reduced_ids)]
  #rownames(rpkm)<-data_rpkm$Hugo_Symbol
  
  #----------------------------------------------------------------------------------------------
  
  #unique names from reduced 
  cancer_types = unique(str_util(reduced_ids))
  cancer_freqs <- factor(str_util(reduced_ids), levels=cancer_types)
  
  #cleaning 
  #rm(data_patient, data_sample, data_treatment_auc, path, reduced_auc_ids,reduced_mrna_ids, reduced_sample_ids)
  #import dataset
  
  list(auc=auc,rpkm=rpkm)
  }

