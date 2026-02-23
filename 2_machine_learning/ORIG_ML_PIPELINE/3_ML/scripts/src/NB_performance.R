library(e1071)
library(caret)
library(pROC)

base_dir <- "/home/nikita.arya/ensemble_pipeline"
input_file <- file.path(base_dir, "output_data/MLDB_repro.csv")
model_dir  <- file.path(base_dir, "output_models/NB")

output_dir <- file.path(base_dir, "output_data/NB_validation")
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

metrics_outfile <- file.path(output_dir, "NB_metrics.csv")
pred_outfile    <- file.path(output_dir, "NB_validation.csv")

df <- read.csv(input_file, check.names = FALSE)
predictor_cols <- 1:79

metrics_list <- list()
gene_names <- colnames(df)[80:ncol(df)]
val_preds_matrix <- matrix(NA, nrow = nrow(df), ncol = length(gene_names))
colnames(val_preds_matrix) <- gene_names
rownames(val_preds_matrix) <- rownames(df)

cat(paste("Starting NB Validation for", length(gene_names), "genes...\n"))

for (i in seq_along(gene_names)) {
  
  gene_id <- i 
  target_name <- gene_names[i]
  
  data <- df[, c(predictor_cols, 79 + i)]
  colnames(data)[80] <- "gene"
  
  set.seed(42)
  trainIndex <- createDataPartition(data$gene, p = 0.8, list = FALSE, times = 1)
  
  val_indices <- setdiff(1:nrow(df), trainIndex)
  dataVal <- data[val_indices, ]
  
  dataVal[] <- lapply(dataVal, function(x) factor(as.character(x), levels = c("0", "1")))
  
  model_file <- file.path(model_dir, paste0("NB", gene_id, ".RData"))
  
  if (!file.exists(model_file)) {
    next
  }
  
  load(model_file) 
  
  tryCatch({
    pred_raw <- predict(model, dataVal[, predictor_cols], type = "raw")
    
    if ("1" %in% colnames(pred_raw)) {
      probs <- pred_raw[, "1"]
    } else {
      probs <- rep(0, nrow(dataVal))
    }
    
    val_preds_matrix[val_indices, i] <- probs
    
    actual_labels <- dataVal$gene
    
    roc_val <- NA
    actual_num <- as.numeric(as.character(actual_labels))
    if (var(actual_num) > 0) {
       roc_obj <- roc(actual_num, probs, quiet = TRUE)
       roc_val <- as.numeric(auc(roc_obj))
    }
    
    pred_class <- factor(ifelse(probs > 0.5, "1", "0"), levels = c("0", "1"))
    cm <- confusionMatrix(pred_class, actual_labels, positive = "1")
    
    metrics_list[[length(metrics_list) + 1]] <- data.frame(
      Gene_ID = gene_id,
      Gene_Name = target_name,
      AUC = roc_val,
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

cat("NB Validation Complete.\n")
