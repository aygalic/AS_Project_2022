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

#ORGANIZES THE MAIN DATASET INTO BLOCKS BASED ON CANCER TYPES 

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

#---------------------------------------------------------------------------------------------
#MUTLICLASS CANCER CELL PCA, SAMPLES ARE DRUGS

querried_names = c("BONE", "SKIN")
sub_auc<- block_dat(querried_names, reduced_auc)
dim_<-dim(sub_auc)

na_fill = 0.5
sub_auc[is.na(sub_auc)]= na_fill
#sum(is.na(sub_auc))
pc <- princomp(sub_auc,scores=T)
#summary(pc)

#IDEA: different drugs should perform optimally on different cell lines
#Pca on auc as features might reveal clusters of drugs which perform similarly on the same cancer types
#TODO: examine loadings, add color map based on argmax(mean(auc on celltype)) over celltypes

#classify what 'kind' of drug we have
drugs = data.frame(rowMeans(sub_auc[,grepl(querried_names[1], colnames(sub_auc))]))
colnames(drugs)<-querried_names[1];

for(c in querried_names[-1]){
  averaged_col = data.frame(rowMeans(sub_auc[,grepl(c, colnames(sub_auc))]));
  colnames(averaged_col)<-c;
  drugs<-cbind(drugs,averaged_col);
}

#extract name of cell type - naive 
#add an indeterminate class when distinction not clear (arbitrary threshold)
#WILL THROW ERROR IF WE ONLY CONSIDER 1 CANCER CLASS

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

drug_max = c(NULL)

for(i in 1:dim_[1]){
  drug_max<-cbind(drug_max,which_cell(drugs[i,], thresh=0.05))
}
drug_max<-t(data.frame(drug_max));
colnames(drug_max)<-"DRUG_MAX"
drugs<-cbind(drugs,drug_max)

#get some colors in here 
classes<-append(querried_names, "INDETERMINATE")
col.cells <- rainbow(length(classes))
col.labels <- rep(NA, dim_[1])

for(i in 1:length(col.labels)){
  col.labels[i] = col.cells[match(drugs[i,"DRUG_MAX"],classes)]
}

#PLOT
x11()
plot(pc$scores[,1], pc$scores[,2],pch=16, col=col.labels)

#----------------------------------------------------------------------------------------------
#MULTICLASS CANCER PCA, SAMPLES ARE CANCERS
#ISSUE: NEED MORE DATA SO n>p
querried_names = c("CENTRAL_NERVOUS_SYSTEM", "BREAST", "BONE", "LIVER", "PANCREAS", "STOMACH")
sub_auc2<- block_dat(querried_names, reduced_auc)

na_fill = 0.5
sub_auc2[is.na(sub_auc2)]= na_fill
sub_auc2<-t(sub_auc2);
#sum(is.na(sub_auc))
pc2 <- princomp(sub_auc2,scores=T)

x11()
plot(pc2$scores[,1], pc2$scores[,2],pch=16)
