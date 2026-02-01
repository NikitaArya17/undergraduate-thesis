library(e1071)

df = read.csv("/home/nikita.arya/ensemble_pipeline/output_data/MLDB_repro.csv",check.names = FALSE)
condition =  read.csv("/home/nikita.arya/ensemble_pipeline/output_data/CultureCondition.csv",check.names = FALSE,
                      stringsAsFactors = FALSE)
ss = paste(colnames(condition), condition[1,], sep= "_")

preds =  as.numeric()
for (gene_id in 1:1990) {
  data=df[,c(1:79,79+gene_id)]
  colnames(data)[80] = "gene"
  data[] <- lapply(data, factor)
  data$gene <- factor(data$gene, levels = c("0", "1"))
  
  test = data[1,1:80,drop=FALSE]
  test[1,] = 0
  test[,ss] = 1
  hyper = c(0.77,0.08)
  model = svm(gene~.,data=data,kernel="radial",cost = hyper[1],
                gamma = hyper[2],probability = TRUE, type = "C")
  save(model,file = paste("/home/nikita.arya/ensemble_pipeline/output_models/SVM",gene_id,".RData",sep=""))
  
  pred = attr(predict(model,test,probability = TRUE),"probabilities")[,'1']
  preds[gene_id] = pred
  if (gene_id %% 100 == 0) gc()
}

res = data.frame(colnames(df)[79+1:1990], preds)
colnames(res) = c("genome site", "probability")
write.csv(res, "/home/nikita.arya/ensemble_pipeline/output_predictions/SVM/pred_SVM.csv", row.names = FALSE)
