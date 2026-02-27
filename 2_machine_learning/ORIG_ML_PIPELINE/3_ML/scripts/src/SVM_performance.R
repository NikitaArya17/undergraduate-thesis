library(e1071)
library(caret)
library(pROC)
library(PRROC)

base_dir <- "/home/nikita.arya/ensemble_pipeline"
input_file <- file.path(base_dir, "output_data/MLDB_repro.csv")
model_dir  <- file.path(base_dir, "output_models/SVM")

output_dir <- file.path(base_dir, "output_data/SVM_validation")
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

metrics_outfile <- file.path(output_dir, "SVM_metrics.csv")
pred_outfile    <- file.path(output_dir, "SVM_validation.csv")

df <- read.csv(input_file, check.names = FALSE)
predictor_cols <- 1:79

metrics_list <- list()

gene_names <- colnames(df)[80:ncol(df)]
val_preds_matrix <- matrix(NA, nrow = nrow(df), ncol = length(gene_names))
colnames(val_preds_matrix) <- gene_names
rownames(val_preds_matrix) <- rownames(df)

cat(paste("Starting SVM Validation for", length(gene_names), "genes...\n"))

for (i in seq_along(gene_names)) {
  
  gene_id <- i
  target_name <- gene_names[i]
  
  data <- df[, c(predictor_cols, 79 + i)]
  colnames(data)[80] <- "gene"
  
  set.seed(2021)
  trainIndex <- createDataPartition(data$gene, p = 0.8, list = FALSE, times = 1)
  
  val_indices <- setdiff(1:nrow(df), trainIndex)
  dataVal <- data[val_indices, ]
  
  dataVal$gene <- factor(dataVal$gene, levels = c("0", "1"))
  
  model_file <- file.path(model_dir, paste0("SVM", gene_id, ".RData"))
  
  if (!file.exists(model_file)) {
    next
  }
  
  load(model_file)
  
  tryCatch({
    pred_obj <- predict(model, dataVal[, predictor_cols], probability = TRUE)
    prob_matrix <- attr(pred_obj, "probabilities")
    probs <- prob_matrix[, "1"]
    
    val_preds_matrix[val_indices, i] <- probs
    
    actual_labels <- dataVal$gene
    
    roc_val <- NA
    prc_val <- NA
    if (length(unique(actual_labels)) > 1) {
       roc_obj <- roc(as.numeric(actual_labels) - 1, probs, quiet = TRUE)
       roc_val <- as.numeric(auc(roc_obj))
       prc_obj <- pr.curve(scores.class0 = probs, weights.class0 = as.numeric(actual_labels) - 1, curve = TRUE)
       prc_val <- as.numeric(prc_obj$auc.integral)
    }
    
    pred_class <- factor(ifelse(probs > 0.5, "1", "0"), levels = c("0", "1"))
    cm <- confusionMatrix(pred_class, actual_labels, positive = "1")
    
    metrics_list[[length(metrics_list) + 1]] <- data.frame(
      Gene_ID = gene_id,
      Gene_Name = target_name,
      AUC = roc_val,
      AUPRC = prc_val,
      Accuracy = cm$overall['Accuracy'],
      Sensitivity = cm$byClass['Sensitivity'],
      Specificity = cm$byClass['Specificity']
    )
    
  }, error = function(e) {
    cat(paste("Error validating gene", gene_id, ":", e$message, "\n"))
  })

  if (i %% 50 == 0) {
    cat(paste("Validated:", i, "\n"))
    gc()
  }
}

cat("Saving metrics...\n")
metrics_df <- do.call(rbind, metrics_list)
write.csv(metrics_df, metrics_outfile, row.names = FALSE)

cat("Saving predictions...\n")
write.csv(val_preds_matrix, pred_outfile, row.names = TRUE)

cat("SVM Validation Complete.\n")
