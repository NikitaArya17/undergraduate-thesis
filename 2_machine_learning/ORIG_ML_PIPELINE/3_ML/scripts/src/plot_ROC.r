library(pROC)
library(ggplot2)

base_dir   <- "3_ML/output_data"
truth_file <- file.path(base_dir, "MLDB_repro.csv")
svm_file   <- file.path(base_dir, "SVM_validation/SVM_validation.csv")
nb_file    <- file.path(base_dir, "NB_validation/NB_validation.csv")
ann_file   <- file.path(base_dir, "ANN_validation/ANN_validation.csv")
output_file <- file.path(base_dir, "Ensemble_ROC.png")

truth_all <- read.csv(truth_file, check.names=FALSE)

svm <- read.csv(svm_file, check.names=FALSE)
nb  <- read.csv(nb_file,  check.names=FALSE)
ann <- read.csv(ann_file, check.names=FALSE)

set.seed(42)
train_idx <- sample(1:nrow(truth_all), 0.8 * nrow(truth_all))
val_idx   <- setdiff(1:nrow(truth_all), train_idx)
truth_val <- truth_all[val_idx, ]

target_rows <- nrow(truth_val) 

pad_to_target <- function(df, target_n) {
  current_n <- nrow(df)
  if (current_n < target_n) {
    diff <- target_n - current_n
    cat(sprintf("  - Padding dataframe with %d empty rows (NAs)\n", diff))
    # Create a block of NAs with same columns
    na_block <- df[1:diff, ]
    na_block[,] <- NA
    # Bind
    return(rbind(df, na_block))
  } else if (current_n > target_n) {
    # If a file is somehow larger, we assume it matches the top 'target_n' rows
    # (or you can throw an error)
    cat(sprintf("  - Warning: Truncating dataframe from %d to %d rows\n", current_n, target_n))
    return(df[1:target_n, ])
  }
  return(df)
}

# C. Apply Padding
svm_pad <- pad_to_target(svm, target_rows)
nb_pad  <- pad_to_target(nb,  target_rows)
ann_pad <- pad_to_target(ann, target_rows)

# D. Ensure Column Alignment
# We only use genes present in ALL files (plus Truth)
common_genes <- intersect(colnames(svm), intersect(colnames(nb), colnames(ann)))
common_genes <- intersect(common_genes, colnames(truth_val))
cat(sprintf("Using %d common genes.\n", length(common_genes)))

# Subset to common columns
svm_final <- svm_pad[, common_genes]
nb_final  <- nb_pad[,  common_genes]
ann_final <- ann_pad[, common_genes]
truth_final <- truth_val[, common_genes]

# --- 4. Calculate Ensemble (Partial Averaging) ---
cat("Calculating Ensemble (ignoring NAs)...\n")

# We use a 3D array to average across the 'Model' dimension
# Dimensions: [Rows, Genes, Models]
data_array <- array(
  data = c(as.matrix(svm_final), as.matrix(nb_final), as.matrix(ann_final)),
  dim = c(nrow(svm_final), ncol(svm_final), 3)
)

# Calculate Mean across the 3rd dimension (Models), ignoring NAs
# This handles the case where ANN is NA; the result is Mean(SVM, NB)
ensemble_final <- apply(data_array, c(1, 2), mean, na.rm = TRUE)

# Convert back to data frame
ensemble_final <- as.data.frame(ensemble_final)
colnames(ensemble_final) <- common_genes

# --- 5. Generate ROC Objects ---
cat("Generating ROC objects (Micro-Average)...\n")

get_roc_obj <- function(preds, truth) {
  p_vec <- as.vector(as.matrix(preds))
  t_vec <- as.vector(as.matrix(truth))
  valid <- !is.na(p_vec) & !is.na(t_vec)
  roc(t_vec[valid], p_vec[valid], quiet=TRUE)
}

roc_list <- list(
  "SVM"      = get_roc_obj(svm_final, truth_final),
  "NB"       = get_roc_obj(nb_final,  truth_final),
  "ANN"      = get_roc_obj(ann_final, truth_final),
  "Ensemble" = get_roc_obj(ensemble_final, truth_final)
)

# --- 6. Plotting ---
cat("Plotting...\n")
g <- ggroc(roc_list, size=1) +
  geom_abline(slope=1, intercept=1, linetype="dashed", color="gray50") +
  scale_color_manual(values = c(
    "SVM" = "#3498db", 
    "NB"  = "#2ecc71", 
    "ANN" = "#9b59b6", 
    "Ensemble" = "#e74c3c"
  )) +
  labs(
    title = "Ensemble ROC",
    color = "Model",
    x = "Specificity",
    y = "Sensitivity"
  ) +
  theme_minimal(base_size = 14) +
  theme(
    legend.position = c(0.75, 0.25),
    legend.background = element_rect(fill="white", color="gray90"),
    aspect.ratio = 1
  )

ggsave(output_file, plot=g, width=8, height=8)
print(g)
cat(sprintf("ROC Curve saved to %s\n", output_file))