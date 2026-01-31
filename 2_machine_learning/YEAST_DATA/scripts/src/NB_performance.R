library(pROC)
library(e1071)


base_dir <- "/home/nikita.arya/ensemble_pipeline/YEAST_DATA"
input_dir <- file.path(base_dir, "input_data/MLDB_yeast.csv")
feature_dir <- file.path(base_dir, "output_data/feature_selection/NB")

output_dir <- file.path(base_dir, "output_predictions/NB")
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

pred_output_file <- file.path(output_dir, "NB_validation.csv")
metrics_output_file <- file.path(output_dir, "NB_validation_metrics.csv")

read_features_robust <- function(filepath) {
  if (!file.exists(filepath)) return(NULL)
  raw_lines <- readLines(filepath, warn = FALSE)
  full_text <- paste(raw_lines, collapse = " ")
  clean_text <- gsub(",", " ", full_text) 
  indices <- as.numeric(unlist(strsplit(clean_text, "\\s+")))
  return(indices[!is.na(indices)])
}


cat("Loading Yeast Data...\n")
df <- read.csv(input_dir, row.names = 1, check.names = FALSE)
total_cols <- ncol(df)
pheno_names <- colnames(df)[(total_cols - 7):total_cols]

results_df <- data.frame(row.names = rownames(df))
metrics_list <- list()

for (p_name in pheno_names) {
  results_df[[paste0("Actual_", p_name)]] <- df[[p_name]]
}


for (pheno_id in 0:7) {
  target_name <- pheno_names[pheno_id + 1]
  target_col_idx <- (total_cols - 7) + pheno_id
  
  cat(paste("\n=== Validating Phenotype:", target_name, "===\n"))
  
  for (setting_id in c(47, 48, 49)) {
    
    feature_file <- file.path(feature_dir, paste0("pheno_class_", pheno_id, "_setting_", setting_id))
    surviving_indices <- read_features_robust(feature_file)
    
    if (is.null(surviving_indices) || length(surviving_indices) == 0) {
      cat(paste("  [Skip] Setting", setting_id, ": No features found.\n"))
      next
    }
    
    valid_mask <- !is.na(df[, target_col_idx])
    
    train_data <- df[valid_mask, c(surviving_indices, target_col_idx)]
    
    y_col_name <- colnames(train_data)[ncol(train_data)]
    train_data[[y_col_name]] <- as.factor(train_data[[y_col_name]])
    
    f <- as.formula(paste("`", y_col_name, "` ~ .", sep=""))
    model <- naiveBayes(f, data = train_data, laplace = 0)
    
    test_data_all <- df[, surviving_indices]
    pred_raw_all <- predict(model, test_data_all, type = "raw")
    probs_all <- pred_raw_all[, "1"]
    
    probs_valid <- probs_all[valid_mask]
    actual_valid <- df[valid_mask, target_col_idx]
    
    
    roc_obj <- roc(actual_valid, probs_valid, quiet = TRUE)
    auc_val <- as.numeric(auc(roc_obj))
    
    pred_class <- ifelse(probs_valid > 0.5, 1, 0)
    acc_val <- sum(pred_class == actual_valid) / length(actual_valid)
    
    cat(paste("  -> Set", setting_id, "| Genes:", length(surviving_indices), 
              "| AUC:", round(auc_val, 3), "| Acc:", round(acc_val, 3), "\n"))
    
    col_name <- paste0("Pred_", target_name, "_Set_", setting_id)
    results_df[[col_name]] <- probs_all
    
    metrics_list[[length(metrics_list) + 1]] <- data.frame(
      Phenotype = target_name,
      Setting = setting_id,
      Num_Genes = length(surviving_indices),
      AUC = auc_val,
      Accuracy = acc_val
    )
  }
}

cat("\nSaving results...\n")

write.csv(results_df, pred_output_file, row.names = TRUE)

metrics_df <- do.call(rbind, metrics_list)
write.csv(metrics_df, metrics_output_file, row.names = FALSE)

cat(paste("Done! Files saved:\n1.", pred_output_file, "\n2.", metrics_output_file, "\n"))
