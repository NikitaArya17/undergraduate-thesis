library(e1071)
library(smotefamily)
library(caret)

df = read.csv("/home/nikita.arya/ensemble_pipeline/output_data/MLDB_repro.csv",check.names = FALSE)
condition =  read.csv("/home/nikita.arya/ensemble_pipeline/output_data/CultureCondition.csv",check.names = FALSE,
                      stringsAsFactors = FALSE)
ss = paste(colnames(condition), condition[1,], sep= "_")

preds =  as.numeric()
for (gene_id in 1:1990) {

  data = df[, c(1:79, 79 + gene_id)]
  colnames(data)[80] = "gene"
  
  set.seed(42)
  trainIndex <- createDataPartition(data$gene, p = .8, list = FALSE, times = 1)
  dataTrain <- data[ trainIndex,]
  
  # If we only have 0s, or only 1s, we can't do SMOTE or SVM
  if (length(unique(dataTrain$gene)) < 2) {
    # If only one class is present... 
    preds[gene_id] = 0 
    
    next 
  }
  
  X_train <- dataTrain[, 1:79]
  y_train <- dataTrain[, 80]
  n_pos <- sum(y_train == "1")
  
  if (n_pos >= 2) {
     k_dynamic <- min(5, n_pos - 1)
     # To address class imbalance
     smote_obj <- SMOTE(X = X_train, target = y_train, K = k_dynamic, dup_size = 0)
     dataTrain_balanced <- smote_obj$data
     colnames(dataTrain_balanced)[ncol(dataTrain_balanced)] <- "gene"
  } else {
     dataTrain_balanced <- dataTrain
  }

  dataTrain_balanced$gene <- factor(dataTrain_balanced$gene, levels = c("0", "1"))
  
  if (length(unique(dataTrain_balanced$gene)) < 2) {
      preds[gene_id] = 0
      next
  }

  hyper = c(0.77, 0.08)
  
  possible_error <- tryCatch({
      model = svm(gene ~ ., data = dataTrain_balanced, kernel = "radial", 
                  cost = hyper[1], gamma = hyper[2], probability = TRUE, type = "C-classification")
      
      
      save(model, file = paste("/home/nikita.arya/ensemble_pipeline/output_models/SVM/SVM", gene_id, ".RData", sep = ""))
      
      test_vector = data[1, 1:79, drop = FALSE] 
      test_vector[1, ] = 0       
      valid_ss <- intersect(ss, colnames(test_vector))
      test_vector[, valid_ss] = 1
      
      pred_prob = attr(predict(model, test_vector, probability = TRUE), "probabilities")[, '1']
      preds[gene_id] = pred_prob
      
  }, error = function(e) {
      print(paste("Error in gene", gene_id, ":", e$message))
      preds[gene_id] <<- 0 
  })
  
  if (gene_id %% 100 == 0) {
    print(paste("Processed gene:", gene_id))
    gc()
  }
}

res = data.frame(colnames(df)[79+1:1990], preds)
colnames(res) = c("genome_site", "probability")
write.csv(res, "/home/nikita.arya/ensemble_pipeline/output_predictions/SVM/pred_SVM.csv", row.names = FALSE)
