cat(getwd())
library(e1071)
library(ROCR)
source("/home/nikita.arya/ensemble_pipeline/src/helper/libs/compute_AUC.R")
source("/home/nikita.arya/ensemble_pipeline/src/helper/backward_wrapper.R")
source("/home/nikita.arya/ensemble_pipeline/src/helper/libs/cross_validation.R")

args = commandArgs(TRUE)
gene_id = as.integer(args[1])

index = as.integer(args[2])
settings = read.csv("/home/nikita.arya/ensemble_pipeline/output_data/SVM_NB_setting",stringsAsFactors = FALSE,header=FALSE)
ML_method = settings[index,1]

cat("setting: ",as.matrix(settings[index,]))

hypers = c(settings[index,2], settings[index,3])

# gene_id = 1
# ML_method = "SVM"
# hypers = c(3,50)

df = read.csv("/home/nikita.arya/ensemble_pipeline/output_data/MLDB_repro.csv")
data=df[,c(1:84,84+gene_id)]
colnames(data)[85] = "gene"
data[] <- lapply(data, factor)
res = backward_wrapper(data,ML_method,hypers)
               
out = data.frame(res[[1]], res[[2]])
colnames(out) = c("AUC", "Excluded_index")

out_file = paste("auc_exclude_index_",gene_id,"_",ML_method,"_",index,"_.csv", sep="")

head(out)
cat(file.path("/home/nikita.arya/ensemble_pipeline/output_data/feature_selection",ML_method,out_file))
write.csv(out,file = file.path("/home/nikita.arya/ensemble_pipeline/output_data/feature_selection",ML_method,out_file),row.names=FALSE)
