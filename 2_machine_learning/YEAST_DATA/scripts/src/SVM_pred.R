library(e1071)

base_dir <- "/home/nikita.arya/ensemble_pipeline/YEAST_DATA"
input_dir <- file.path(base_dir, "input_data/MLDB_yeast.csv")
feature_dir <- file.path(base_dir, "output_data/feature_selection/SVM")
pred_dir <- file.path(base_dir, "output_predictions/SVM")
dir.create(pred_dir, recursive = TRUE, showWarnings = FALSE)


get_svm_settings <- function(id) {
  if (id == 47) return(list(cost = 1,   gamma = 0.01))
  if (id == 48) return(list(cost = 10,  gamma = 0.001))
  if (id == 49) return(list(cost = 100, gamma = 0.0001))
  return(list(cost = 1, gamma = 0.01))
}

df <- read.csv(input_dir, row.names = 1, check.names = FALSE)
total_cols <- ncol(df)
pheno_names <- colnames(df)[(total_cols - 7):total_cols]
results_df <- data.frame(row.names = rownames(df))


for (pheno_id in 0:7) {
  target_name <- pheno_names[pheno_id + 1]
  target_col_idx <- (total_cols - 7) + pheno_id

  for (setting_id in c(47, 48, 49)) {
    cat(paste("\n--- Predicting Phenotype:", target_name, "| Setting:", setting_id, "---\n"))

    
    feature_file <- file.path(feature_dir, paste0("pheno_class_", pheno_id, "_setting_", setting_id))

    if (!file.exists(feature_file)) {
      cat("Warning: Feature file missing. Skipping.\n")
      next
    }

    
    surviving_indices <- as.numeric(read.table(feature_file)$V1)

    
    params <- get_svm_settings(setting_id)

    train_mask <- !is.na(df[, target_col_idx])
    train_data <- df[train_mask, c(surviving_indices, target_col_idx)]
    
    
    y_col_name <- colnames(train_data)[ncol(train_data)]
    train_data[[y_col_name]] <- as.factor(train_data[[y_col_name]])

    
    model <- svm(as.formula(paste("`", y_col_name, "` ~ .", sep="")), 
                 data = train_data, 
                 kernel = "radial",
                 cost = params$cost, 
                 gamma = params$gamma,
                 probability = TRUE)

    
    test_data <- df[, surviving_indices]
    pred_obj <- predict(model, test_data, probability = TRUE)
    
    probs <- attr(pred_obj, "probabilities")[, "1"]

    results_df[, paste0("Pred_", target_name, "_Set_", setting_id)] <- probs
  }
  
  
  results_df[, paste0("Actual_", target_name)] <- df[, (total_cols - 7) + pheno_id]
}

write.csv(results_df, file.path(pred_dir, "SVM_pred.csv"), row.names = TRUE)
cat("\nAll SVM predictions saved to SVM_pred.csv\n")
