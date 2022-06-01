#first get rid of nrows line in rpkm dataset import in exp_base_script.R
#also make sure to provide the right path to the dataset

source("../utils/nate_utils.R")
source("../utils/exp_base_script.R")

View(rpkm)

#extract whatever
cancer_types
cancers = c("CENTRAL_NERVOUS_SYSTEM", "BREAST")
data.rpkm = block_dat(cancers, rpkm)

#auc if ya want 
data.auc = block_dat(cancers, auc)

colnames(data.rpkm) == colnames(data.auc)
