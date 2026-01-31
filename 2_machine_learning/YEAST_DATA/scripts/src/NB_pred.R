library(e1071)

base_dir <- "/home/nikita.arya/ensemble_pipeline/YEAST_DATA"
input_dir <- file.path(base_dir, "input_data/MLDB_yeast.csv")
feature_dir <- file.path(base_dir, "output_data/feature_selection/NB")
pred_dir <- file.path(base_dir, "output_predictions/NB")

dir.create(pred_dir, recursive = TRUE, showWarnings = FALSE)

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

for (pheno_id in 0:7) {
  target_name <- pheno_names[pheno_id + 1]
  target_col_idx <- (total_cols - 7) + pheno_id
  
  cat(paste("\n--- Predicting Phenotype:", target_name, "---\n"))
  
  for (setting_id in c(47, 48, 49)) {
    
    feature_file <- file.path(feature_dir, paste0("pheno_class_", pheno_id, "_setting_", setting_id))
    surviving_indices <- read_features_robust(feature_file)
    
    if (is.null(surviving_indices) || length(surviving_indices) == 0) {
      cat(paste("  [Skip] Setting", setting_id, ": No features found.\n"))
      next
    }

    train_mask <- !is.na(df[, target_col_idx])
    train_data <- df[train_mask, c(surviving_indices, target_col_idx)]
    
    y_col_name <- colnames(train_data)[ncol(train_data)]
    train_data[[y_col_name]] <- as.factor(train_data[[y_col_name]])
    
    f <- as.formula(paste("`", y_col_name, "` ~ .", sep=""))
    model <- naiveBayes(f, data = train_data, laplace = 0)
    
    test_data <- df[, surviving_indices]
    
    pred_raw <- predict(model, test_data, type = "raw")
    
    probs <- pred_raw[, "1"]
    
    col_name <- paste0("Pred_", target_name, "_Set_", setting_id)
    results_df[[col_name]] <- probs
    
    cat(paste("  -> Setting", setting_id, ": Prediction complete using", length(surviving_indices), "genes.\n"))
  }
  
  results_df[[paste0("Actual_", target_name)]] <- df[[target_name]]
}

output_file <- file.path(pred_dir, "NB_pred.csv")
write.csv(results_df, output_file, row.names = TRUE)
cat(paste("\nDone! All NB predictions saved to:", output_file, "\n"))
