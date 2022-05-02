#you can include the pre-processing i've done and the utility functions like this 
your_path_to_repo<- "/Users/nathanielnethercott/Desktop/School/Polimi/2021:2022/AS/AS_Project_2022"

source(paste(your_path_to_repo,"/nate/utils/nate_utils.R", sep=""))
source(paste(your_path_to_repo,"/nate/utils/exp_base_script.R", sep=""))

#I still need to cleanup all the environment variables that get imported but the 
#ones you want are 
View(auc)
dim(auc)
View(rpkm)

#I also have a factor of the cancer types on the basis of patient_id 
table(cancer_freqs)

#NOTE: THESE ARE DIFFERENT LABELS THAN FROM THE "Cancer Type Detailed" which we used before

#To access sub dataframe of specific cancers you just use the function block_dat
#e.g.
cancers<-c("URINARY_TRACT")

sub_auc<-block_dat(cancers, auc)
View(sub_auc)
dim(sub_auc)

#with rpkm its a little different.  All the data is stored the same way as auc, 
#but since there's duplicate genes (which you noticed), i couldn't write the rownames of rpkm as 
#the genes.  Instead they're unlabeled but i have a separate variable for the genes 
sub_rpkm <-block_dat(cancers, rpkm)
View(sub_rpkm)
dim(sub_rpkm)

hugo
