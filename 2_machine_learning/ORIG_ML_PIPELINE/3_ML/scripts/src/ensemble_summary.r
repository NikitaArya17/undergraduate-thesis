library(pROC)
library(PRROC)

base_dir <- "3_ML/output_data"
truth_file <- file.path(base_dir, "MLDB_repro.csv")
svm_file   <- file.path(base_dir, "SVM_validation/SVM_validation.csv")
nb_file    <- file.path(base_dir, "NB_validation/NB_validation.csv")
ann_file   <- file.path(base_dir, "ANN_validation/ANN_validation.csv")

truth_all <- read.csv(truth_file, check.names = FALSE)
set.seed(42)
train_idx <- sample(1:nrow(truth_all), 0.8 * nrow(truth_all))
val_idx   <- setdiff(1:nrow(truth_all), train_idx)
truth_val <- truth_all[val_idx, ]

svm_preds <- read.csv(svm_file, check.names = FALSE)
nb_preds  <- read.csv(nb_file,  check.names = FALSE)
ann_preds <- read.csv(ann_file, check.names = FALSE)

ensemble_preds <- (svm_preds + nb_preds + ann_preds) / 3

calc_metrics <- function(preds_df, truth_df, model_name) {
  cat(sprintf("Processing %s...\n", model_name))
  aucs <- c(); auprcs <- c()
  genes <- colnames(preds_df)
  
  for (gene in genes) {
    yp <- preds_df[[gene]]
    yt <- truth_df[[gene]]
    
    valid <- !is.na(yp) & !is.na(yt)
    if (sum(valid) < 5 || length(unique(yt[valid])) < 2) next
    
    yp_v <- yp[valid]; yt_v <- yt[valid]
    
    aucs <- c(aucs, as.numeric(auc(roc(yt_v, yp_v, quiet=TRUE))))
    fg <- yp_v[yt_v == 1]; bg <- yp_v[yt_v == 0]
    pr <- pr.curve(scores.class0 = fg, scores.class1 = bg, curve = FALSE)
    auprcs <- c(auprcs, pr$auc.integral)
  }
  
  fmt <- function(vals) {
    m <- mean(vals, na.rm=TRUE)
    se <- sd(vals, na.rm=TRUE) / sqrt(length(na.omit(vals)))
    return(sprintf("%.4f ± %.4f", m, se))
  }
  
  return(data.frame(
    Model = model_name,
    AUC = fmt(aucs),
    AUPRC = fmt(auprcs),
    "AUC/AUPRC" = sprintf("%.2f / %.2f", mean(aucs), mean(auprcs))
  ))
}

final_results <- rbind(
  calc_metrics(svm_preds, truth_val, "SVM"),
  calc_metrics(nb_preds,  truth_val, "NB"),
  calc_metrics(ann_preds, truth_val, "ANN"),
  calc_metrics(ensemble_preds, truth_val, "Ensemble")
)

print(final_results)
write.csv(final_results, "Ensemble_model_comparison.csv", row.names=FALSE)