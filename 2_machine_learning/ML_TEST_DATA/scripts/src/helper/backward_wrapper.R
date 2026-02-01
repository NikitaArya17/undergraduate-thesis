backward_wrapper = function(data, ML_method, hypers) {
  df = data
  # Target is the last column
  target_col_idx = ncol(df)
  
  # survive_index starts with all genes (all columns except the last one)
  survive_index = seq(1, target_col_idx - 1, 1)
  
  # survive_index = sample(survive_index, 55) 
  
  best_perf = 0
  best_perfs = as.numeric()
  exclude_index = as.numeric()
  
  # Loop until we only have one gene left or performance stops improving
  for (i in 1:(length(survive_index) - 1)) {
    cat("\nIteration:", i, "| Remaining genes:", length(survive_index), "\n")
    perfs = as.numeric()
    
    for(idx in survive_index) {
      test_indices = survive_index[survive_index != idx]
      
      # We create a temporary dataframe with only the surviving genes + target
      temp_df = df[, c(test_indices, target_col_idx)]
      
      # Run cross validation
      res = cross_validation(temp_df, ML_method, hypers)
      
      # Ensure the helper function actually returns "perf"
      
      val = if(!is.null(res[["perf"]])) res[["perf"]] else res[[1]]
      perfs = c(perfs, val)
    }
    
    # Identify which gene's removal resulted in the best performance
    max_perf_idx = which.max(perfs)
    current_best = perfs[max_perf_idx]
    
    # If removing the gene helped (or didn't hurt much), exclude it
    if (current_best >= best_perf && best_perf < 0.98) {
      best_perf = current_best
      best_perfs = c(best_perfs, best_perf)
      
      # Identify the actual gene index to drop
      dropped_gene_idx = survive_index[max_perf_idx]
      exclude_index = c(exclude_index, dropped_gene_idx)
      
      # Update survive_index by removing the dropped gene
      survive_index = survive_index[-max_perf_idx]
      
      cat("Best AUC so far:", best_perf, "| Dropped gene index:", dropped_gene_idx, "\n")
    } else {
      cat("Performance plateaued. Stopping selection.\n")
      break
    }
  }
  
  return(list(best_perfs, exclude_index))
}
