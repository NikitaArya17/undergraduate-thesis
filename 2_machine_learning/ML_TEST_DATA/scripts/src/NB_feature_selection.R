library(e1071)

args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 2) {
  stop("Error: Missing arguments. Usage: Rscript NB_feat_select.R <pheno_idx> <setting_id>")
}

pheno_idx <- as.numeric(args[1])
setting_id <- as.numeric(args[2])

get_nb_settings <- function(id) {
  return(list(laplace = 0)) 
}

params <- get_nb_settings(setting_id)

base_dir <- "/home/nikita.arya/ensemble_pipeline/YEAST_DATA"
input_dir <- file.path(base_dir, "input_data/MLDB_yeast.csv")
output_dir <- file.path(base_dir, "output_data/feature_selection/NB")
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

df <- read.csv(input_dir, row.names = 1)

n_phenos <- 8
target_col_idx <- (ncol(df) - n_phenos + 1) + pheno_idx
target_name <- colnames(df)[target_col_idx]

valid_rows <- !is.na(df[[target_name]])
train_data <- df[valid_rows, ]
y <- as.factor(train_data[[target_name]])

current_genes <- 1:100

cat(paste("Starting NB Selection:", target_name, "| Setting:", setting_id, "\n"))

while(length(current_genes) > 10) {

  scores <- sapply(1:length(current_genes), function(i) {
    test_indices <- current_genes[-i]

    
    fit <- naiveBayes(train_data[, test_indices], y, laplace = params$laplace)

    pred <- predict(fit, train_data[, test_indices], type = "class")
    
    return(sum(pred == y) / length(y))
  })

  best_to_remove <- which.max(scores)
  current_genes <- current_genes[-best_to_remove]

  if (length(current_genes) %% 10 == 0) {
    cat(paste("  -", target_name, "(Set", setting_id, "):", length(current_genes), "genes remaining...\n"))
  }
}

output_file <- file.path(output_dir, paste0("pheno_class_", pheno_idx, "_setting_", setting_id))

write.table(current_genes, file = output_file, row.names = FALSE, col.names = FALSE)

cat(paste("SUCCESS: Saved to", output_file, "\n"))
