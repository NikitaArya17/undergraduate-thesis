base_dir <- "3_ML/SMOTE_NB_SVM_only/random_seed_42"

svm_file      <- file.path(base_dir, "SVM_metrics.csv")
nb_file       <- file.path(base_dir, "NB_metrics.csv")
ann_file      <- file.path(base_dir, "ANN_metrics.csv")
ensemble_file <- file.path(base_dir, "Ensemble_metrics.csv")

comparison_outfile <- file.path(base_dir, "Ensemble_model_comparison.csv")

calc_stats <- function(values) {
  v <- as.numeric(na.omit(values))
  n <- length(v)
  
  if (n < 2) return(c(Mean=mean(v), SD=0, SE=0))
  
  m  <- mean(v)
  s  <- sd(v)
  se <- s / sqrt(n)
  
  return(c(Mean=m, SD=s, SE=se))
}

get_model_stats <- function(filepath, model_name) {
  
  if (!file.exists(filepath)) return(NULL)
  cat(sprintf("Processing %s...\n", model_name))
  
  df <- read.csv(filepath, stringsAsFactors = FALSE)
  
  df <- df[!df[[1]] %in% c("AVERAGE", "STD_DEV"), ]
  
  colnames(df) <- tolower(colnames(df))
  
  target_metrics <- list(
    "AUC" = c("auc", "roc_auc"),
    "AUPRC" = c("auprc", "auprc_val"),
    "Accuracy" = c("accuracy", "acc"),
    "Precision" = c("precision", "prec"),
    "Recall" = c("sensitivity", "sens", "recall"),
    "F1" = c("f1_score", "f1")
  )

  stats_row <- data.frame(Model = model_name, stringsAsFactors = FALSE)
  
  for (metric_name in names(target_metrics)) {
    possible_names <- target_metrics[[metric_name]]
    match_col <- intersect(possible_names, colnames(df))
    
    if (length(match_col) > 0) {
      vals <- df[[match_col[1]]]
      stats <- calc_stats(vals)
      
      stats_row[[paste0(metric_name, "_Mean")]] <- stats["Mean"]
      stats_row[[paste0(metric_name, "_SE")]]   <- stats["SE"]
      
    } else {
      stats_row[[paste0(metric_name, "_Mean")]] <- 0
      stats_row[[paste0(metric_name, "_SE")]]   <- 0
    }
  }
  
  return(stats_row)
}

s_svm <- get_model_stats(svm_file, "SVM")
s_nb  <- get_model_stats(nb_file, "NB")
s_ann <- get_model_stats(ann_file, "ANN")
s_ens <- get_model_stats(ensemble_file, "Ensemble")

final_df <- rbind(s_svm, s_nb, s_ann, s_ens)

num_cols <- sapply(final_df, is.numeric)
final_df[num_cols] <- round(final_df[num_cols], 4)

write.csv(final_df, comparison_outfile, row.names = FALSE)

cat("Done. Results saved.")