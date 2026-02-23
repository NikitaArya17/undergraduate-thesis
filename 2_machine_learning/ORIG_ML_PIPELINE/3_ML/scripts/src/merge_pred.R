ANN_pred = read.csv("3_ML/output_pred/pred_ANN.csv")
SVM_pred = read.csv("3_ML/output_pred/pred_ANN.csv")
NB_pred = read.csv("3_ML/output_pred/pred_ANN.csv")

preds = ANN_pred
preds[,2] = (1/3)*(ANN_pred[,2]+SVM_pred[,2]+NB_pred[,2])
colnames(preds) = c("genome site", "probability")
write.csv(preds, "3_ML/output_pred/pred_average.csv", row.names = FALSE)