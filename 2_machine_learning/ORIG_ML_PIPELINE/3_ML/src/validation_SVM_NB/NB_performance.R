require(PRROC)
library(e1071)
library(ROCR)

source("/home/nikita.arya/ensemble_pipeline/src/helper/libs/compute_AUC_PR.R")
source("/home/nikita.arya/ensemble_pipeline/src/helper/backward_wrapper.R")
source("/home/nikita.arya/ensemble_pipeline/src/helper/libs/cross_validation_raw_res.R")
load("home/nikita.arya/ensemble_pipeline/output_data/feature_selection/selected_features/NB_features.rds")

df = read.csv("home/nikita.arya/ensemble_pipeline/output_data/MLDB.csv")

score = as.numeric()

args <- commandArgs(TRUE)
repeat_id  <- args[1]

gene_num = 1990

for (gene_id in 1:gene_num) {

    data=df[,c(1:84,84+gene_id)]
    colnames(data)[85] = "gene"
    data[] <- lapply(data, factor)
    pp = ncol(data)
    index = sample(seq_len(nrow(data)),nrow(data)*2,replace = TRUE)
    train = data[index,]
    test = data[2,-pp,drop=FALSE]
    model = naiveBayes(gene ~ ., train)
    score[gene_id] = predict(model,test,type='raw')[,2]
    
}

out = data.frame(colnames(df)[85:(84+gene_num)], score)
colnames(out) = c("gene","pred")
write.csv(out, file.path("home/nikita.arya/ensemble_pipeline/output_data","NB_validation",paste("prediction_",repeat_id,".csv",sep="")), row.names = FALSE)
