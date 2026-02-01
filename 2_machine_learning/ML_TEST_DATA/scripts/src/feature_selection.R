cat(getwd())
library(e1071)
library(ROCR)

base_dir <- "/home/nikita.arya/ensemble_pipeline/YEAST_DATA"
source(file.path(base_dir, "src/helper/libs/compute_AUC.R"))
source(file.path(base_dir, "src/helper/backward_wrapper.R"))
source(file.path(base_dir, "src/helper/libs/cross_validation.R"))

args = commandArgs(TRUE)
pheno_id = as.integer(args[1])

index = as.integer(args[2])
settings = read.csv(file.path(base_dir, "output_data/SVM_NB_setting"),stringsAsFactors = FALSE,header=FALSE)
ML_method = settings[index,1]

cat("setting: ",as.matrix(settings[index,]))

hypers = c(settings[index,2], settings[index,3])

# gene_id = 1
# ML_method = "SVM"
# hypers = c(3,50)

df = read.csv(file.path(base_dir, "input_data/MLDB_yeast.csv"), row.names = 1, check.names = FALSE)
total_cols <- ncol(df)
gene_cols <- 1:(total_cols - 8)
target_col_id <- (total_cols - 8) + pheno_id

data_subset <- df[, c(gene_cols, target_col_id)]
colnames(data_subset)[ncol(data_subset)] <- "pheno_class"

data_subset <- data_subset[!is.na(data_subset$pheno_class), ]

data_subset$pheno_class <- as.factor(data_subset$pheno_class)

variances <- apply(data_subset[, 1:length(gene_cols)], 2, var)
valid_genes <- which(variances > 0)

cat(paste("Screening", length(valid_genes), "variable genes...\n"))

p_values <- sapply(valid_genes, function(i) {
  tryCatch({
    t.test(data_subset[,i] ~ data_subset$pheno_class)$p.value
  }, error = function(e) return(1.0)) # Assign p=1 to problematic genes
})

top_100_indices <- valid_genes[order(p_values)[1:100]]

final_data <- data_subset[, c(top_100_indices, ncol(data_subset))]
res = backward_wrapper(final_data,ML_method,hypers)
               
out = data.frame(res[[1]], res[[2]])
colnames(out) = c("AUC", "Excluded_index")

out_file = paste("auc_exclude_index_pheno", pheno_id, "_", ML_method, "_", index, ".csv", sep="")
write.csv(out, file = file.path(base_dir, "output_data/feature_selection", ML_method, out_file), row.names=FALSE)

cat("Done successfully for", ML_method, "phenotype", pheno_id)
