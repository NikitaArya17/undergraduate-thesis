library(e1071)
library(smotefamily)
library(caret)

df = read.csv("/home/nikita.arya/ensemble_pipeline/output_data/MLDB_repro.csv", check.names = FALSE)
condition = read.csv("/home/nikita.arya/ensemble_pipeline/output_data/CultureCondition.csv", check.names = FALSE, stringsAsFactors = FALSE)

ss = paste(colnames(condition), condition[1,], sep= "_")

preds = as.numeric()

dir.create("/home/nikita.arya/ensemble_pipeline/output_models/NB", recursive = TRUE, showWarnings = FALSE)
dir.create("/home/nikita.arya/ensemble_pipeline/output_predictions/NB", recursive = TRUE, showWarnings = FALSE)

for (gene_id in 1:1990) {
  
  data = df[, c(1:79, 79 + gene_id)]
  colnames(data)[80] = "gene"
  
  set.seed(42)
  trainIndex <- createDataPartition(data$gene, p = 0.8, list = FALSE, times = 1)
  dataTrain <- data[ trainIndex,]
  
  if (length(unique(dataTrain$gene)) < 2) {
    preds[gene_id] = 0 
    next 
  }
  
  X_train <- dataTrain[, 1:79]
  y_train <- dataTrain[, 80]
  n_pos <- sum(y_train == 1)
  
  if (n_pos >= 2) {
     k_dynamic <- min(5, n_pos - 1)
     smote_obj <- SMOTE(X = X_train, target = y_train, K = k_dynamic, dup_size = 0)
     dataTrain_balanced <- smote_obj$data
     colnames(dataTrain_balanced)[ncol(dataTrain_balanced)] <- "gene"
  } else {
     dataTrain_balanced <- dataTrain
  }
  
  if (length(unique(dataTrain_balanced$gene)) < 2) {
      preds[gene_id] = 0
      next
  }
  numeric_cols <- sapply(dataTrain_balanced, is.numeric)
    dataTrain_balanced[numeric_cols] <- lapply(dataTrain_balanced[numeric_cols], round)
    
    dataTrain_balanced[] <- lapply(dataTrain_balanced, function(x) factor(as.character(x), levels = c("0", "1")))
    
    if (length(unique(dataTrain_balanced$gene)) < 2) {
        preds[gene_id] = 0
        next
    }
  tryCatch({
      model = naiveBayes(gene ~ ., data = dataTrain_balanced)
      save(model, file = paste("/home/nikita.arya/ensemble_pipeline/output_models/NB/NB", gene_id, ".RData", sep = ""))
      
      test_vector = data[1, 1:80, drop = FALSE] 
      test_vector[1, ] = 0       
      valid_ss <- intersect(ss, colnames(test_vector))
      test_vector[, valid_ss] = 1
      
      test_vector[] <- lapply(test_vector, function(x) factor(x, levels = c("0", "1")))

      pred_prob = predict(model, test_vector[, -80, drop=FALSE], type = 'raw')[, "1"]
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

res = data.frame(colnames(df)[79 + 1:1990], preds)
colnames(res) = c("genome_site", "probability")
write.csv(res, "/home/nikita.arya/ensemble_pipeline/output_predictions/NB/pred_NB.csv", row.names = FALSE)
