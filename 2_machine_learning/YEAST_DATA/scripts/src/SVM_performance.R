library(e1071)
library(pROC)

base_dir <- "/home/nikita.arya/ensemble_pipeline/YEAST_DATA"
input_dir <- file.path(base_dir, "input_data/MLDB_yeast.csv")
feature_dir <- file.path(base_dir, "output_data/feature_selection/SVM")

pred_dir <- file.path(base_dir, "output_predictions/SVM")
dir.create(pred_dir, recursive = TRUE, showWarnings = FALSE)

pred_output_file <- file.path(pred_dir, "SVM_validation.csv")
metrics_output_file <- file.path(pred_dir, "SVM_metrics_summary.csv")

read_features_robust <- function(filepath) {
  raw_lines <- readLines(filepath, warn = FALSE)
  full_text <- paste(raw_lines, collapse = " ")
  clean_text <- gsub(",", " ", full_text)
  indices <- as.numeric(unlist(strsplit(clean_text, "\\s+")))
  indices <- indices[!is.na(indices)]
  return(indices)
}

get_svm_settings <- function(id) {
  if (id == 47) return(list(cost = 1,   gamma = 0.01))
  if (id == 48) return(list(cost = 10,  gamma = 0.001))
  if (id == 49) return(list(cost = 100, gamma = 0.0001))
  return(list(cost = 1, gamma = 0.01)) # Default
}

cat("Loading Yeast Data...\n")
df <- read.csv(input_dir, row.names = 1, check.names = FALSE)
total_cols <- ncol(df)

pheno_names <- colnames(df)[(total_cols - 7):total_cols]

results_df <- data.frame(row.names = rownames(df))
metrics_list <- list() # To store rows of metrics


for (p_name in pheno_names) {
  results_df[[paste0("Actual_", p_name)]] <- df[[p_name]]
}

for (pheno_id in 0:7) {
  target_name <- pheno_names[pheno_id + 1]
  target_col_idx <- (total_cols - 7) + pheno_id
  
  cat(paste("\n=== Processing Phenotype:", target_name, "===\n"))
  
  for (setting_id in c(47, 48, 49)) {
    
    feature_file <- file.path(feature_dir, paste0("pheno_class_", pheno_id, "_setting_", setting_id))
    
    if (!file.exists(feature_file)) {
      cat(paste("  [Skip] Missing file for Setting", setting_id, "\n"))
      next
    }
    
    surviving_indices <- read_features_robust(feature_file)
    
    if (length(surviving_indices) == 0) {
      cat(paste("  [Skip] No genes found in file for Setting", setting_id, "\n"))
      next
    }
    
    train_mask <- !is.na(df[, target_col_idx])
    train_data <- df[train_mask, c(surviving_indices, target_col_idx)]
    
    y_col_name <- colnames(train_data)[ncol(train_data)]
    train_data[[y_col_name]] <- as.factor(train_data[[y_col_name]])
    
    params <- get_svm_settings(setting_id)
    
    f <- as.formula(paste("`", y_col_name, "` ~ .", sep=""))
    
    model <- svm(f, data = train_data, 
                 kernel = "radial", 
                 cost = params$cost, 
                 gamma = params$gamma, 
                 probability = TRUE)
    
    test_data <- df[, surviving_indices]
    pred_obj <- predict(model, test_data, probability = TRUE)
    
    probs <- attr(pred_obj, "probabilities")[, "1"]
    
    col_name <- paste0("Pred_", target_name, "_Set_", setting_id)
    results_df[[col_name]] <- probs
    
    
    val_probs <- probs[train_mask]
    val_actual <- df[train_mask, target_col_idx]
    
    
    roc_obj <- roc(val_actual, val_probs, quiet = TRUE)
    auc_val <- as.numeric(auc(roc_obj))
    
    pred_class <- ifelse(val_probs > 0.5, 1, 0)
    acc_val <- sum(pred_class == val_actual) / length(val_actual)
    
    cat(paste("  -> Set", setting_id, "| Genes:", length(surviving_indices), 
              "| AUC:", round(auc_val, 3), "| Acc:", round(acc_val, 3), "\n"))
    
    metrics_list[[length(metrics_list) + 1]] <- data.frame(
      Phenotype = target_name,
      Setting = setting_id,
      Num_Genes = length(surviving_indices),
      AUC = auc_val,
      Accuracy = acc_val
    )
  }
}

cat("\nSaving predictions...\n")
write.csv(results_df, pred_output_file, row.names = TRUE)

cat("Saving metrics summary...\n")
metrics_df <- do.call(rbind, metrics_list)
write.csv(metrics_df, metrics_output_file, row.names = FALSE)

cat(paste("Done! Results saved to:\n1.", pred_output_file, "\n2.", metrics_output_file, "\n"))
