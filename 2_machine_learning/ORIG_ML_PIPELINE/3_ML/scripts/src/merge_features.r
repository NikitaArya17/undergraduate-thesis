library(data.table)

args <- commandArgs(trailingOnly = TRUE)
index <- 0

if (length(args) == 0) {
  stop("Please provide an ML method.")
} else if (!(args[1] %in% c("NB", "SVM", "ANN"))) {
  stop("Invalid method. Please select one of: NB, SVM, ANN")
} else {
  ML_method <- args[1]
}

index <- if (ML_method == "NB") 1 else 2

work_dir <- file.path("/home/nikita.arya/ensemble_pipeline/output_data/feature_selection", ML_method)
output_path <- file.path(work_dir, paste0(ML_method, "_perf_feature_selection.csv"))

if (ML_method != "ANN") {
  file_list <- paste0(work_dir, "/auc_exclude_index_", 1:1990, "_", ML_method, "_", index, "_", ".csv")
  valid_files <- file_list[file.exists(file_list)]
  data_list <- lapply(valid_files, read.csv)
  final_data <- do.call(rbind, data_list)
  write.csv(final_data, output_path, row.names = FALSE)
} else {
  combos <- expand.grid(gene = 1:1990, setting = 47:49)
  file_list <- paste0(work_dir, "/gene_", combos$gene, "_setting_", combos$setting)
  valid_files <- file_list[file.exists(file_list)]
  data_list <- lapply(valid_files, function(f) { if (file.info(f)$size > 0) { fread(f, header = FALSE, fill = TRUE) }})
  final_data <- rbindlist(data_list, fill = TRUE)
  fwrite(final_data, output_path)
}

