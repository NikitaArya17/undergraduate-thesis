library(pROC)
library(MLmetrics)

base_dir <- "3_ML/SMOTE_NB_SVM_only/random_seed_123"

svm_file <- file.path(base_dir, "SVM_validation.csv")
nb_file <- file.path(base_dir, "NB_validation.csv")
ann_file <- file.path(base_dir, "ANN_validation.csv")
truth_file <- file.path(base_dir, "MLDB_repro.csv")

metrics_outfile <- file.path(base_dir, "Ensemble_metrics.csv")
preds_outfile <- file.path(base_dir, "pred_Ensemble.csv")

svm <- read.csv(svm_file, na.strings = c("", "NA"), row.names = 1)
nb <- read.csv(nb_file, na.strings = c("", "NA"), row.names = 1)
ann <- read.csv(ann_file, na.strings = c("", "NA"), row.names = 1)

dim(svm)
dim(ann)
dim(nb)

ensemble_preds <- (svm + nb + ann) / 3

full_data <- read.csv(truth_file, check.names = FALSE)

results <- data.frame(
  Gene = character(),
  Accuracy = numeric(),
  AUC = numeric(),
  AUPRC = numeric(),
  Sensitivity = numeric(),
  Specificity = numeric(),
  Precision = numeric(),
  F1_Score = numeric(),
  MCC = numeric(),
  stringsAsFactors = FALSE
)

genes <- colnames(ensemble_preds)
total_genes <- length(genes)

for (i in 1:total_genes) {
  gene_name <- genes[i]
  
  preds_all <- ensemble_preds[, gene_name]
  
  valid_idx <- which(!is.na(preds_all))
  
  if (length(valid_idx) == 0) next
  
  y_prob <- preds_all[valid_idx]
  y_class <- ifelse(y_prob > 0.5, 1, 0)
  
  if (gene_name %in% colnames(full_data)) {
    y_true <- full_data[valid_idx, gene_name]
  } else {
    target_idx <- 79 + i 
    if (target_idx <= ncol(full_data)) {
        y_true <- full_data[valid_idx, target_idx]
    } else {
        next
    }
  }
  
  cm <- table(factor(y_class, levels=c(0,1)), factor(y_true, levels=c(0,1)))
  tn <- cm[1,1]
  fp <- cm[2,1]
  fn <- cm[1,2]
  tp <- cm[2,2]
  
  acc <- (tp + tn) / sum(cm)
  
  auc_val <- 0.5
  auprc_val <- 0.5

  tryCatch({
    if (length(unique(y_true)) > 1) {
      roc_obj <- roc(y_true, y_prob, quiet=TRUE)
      auc_val <- as.numeric(auc(roc_obj))
    }
  }, error = function(e) {auc_val <<- 0.5 })

  tryCatch({
    if (length(unique(y_true)) > 1) {
      auprc_val <- PRAUC(y_prob, y_true)
    }
  }, error = function(e) {auprc_val <<- 0.5 })

  sens <- ifelse((tp + fn) > 0, tp / (tp + fn), 0)
  spec <- ifelse((tn + fp) > 0, tn / (tn + fp), 0)
  prec <- ifelse((tp + fp) > 0, tp / (tp + fp), 0)
  f1 <- ifelse((prec + sens) > 0, 2 * (prec * sens) / (prec + sens), 0)

  numerator <- (tp * tn) - (fp * fn)
  denominator <- sqrt(as.numeric(tp + fp) * as.numeric(tp + fn) * as.numeric(tn + fp) * as.numeric(tn + fn))
  mcc <- ifelse(denominator == 0, 0, numerator / denominator)
  
  new_row <- data.frame(gene_name, acc, auc_val, auprc_val, sens, spec, prec, f1, mcc, stringsAsFactors=FALSE)
  names(new_row) <- names(results)
  results <- rbind(results, new_row)

}

means <- colMeans(results[,-1], na.rm = TRUE)
sds <- apply(results[,-1], 2, sd, na.rm = TRUE)

summary_mean_df <- data.frame(Gene = "Mean", as.list(means), stringsAsFactors=FALSE)
summary_sd_df <- data.frame(Gene = "Std_Dev", as.list(sds),   stringsAsFactors=FALSE)

colnames(summary_mean_df) <- colnames(results)
colnames(summary_sd_df) <- colnames(results)

results_final <- rbind(results, summary_mean_df, summary_sd_df)

write.csv(results_final, metrics_outfile, row.names = FALSE)
write.csv(ensemble_preds, preds_outfile, row.names = FALSE)
cat("Done. Results saved.")
