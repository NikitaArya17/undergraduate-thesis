library(dplyr)
library(readr)

mut_un <- read_csv("1_create_DB/output_data/mutation/mutation_unified.csv")
one_hot <- read_csv("1_create_DB/output_data/mutation/one_hot.csv")
sort_del <- read_csv("1_create_DB/output_data/mutation/sort_delete_mutation.csv")
sort_rep <- read_csv("1_create_DB/output_data/mutation/sort_replicate.csv")
syn_rep <- read_csv("1_create_DB/output_data/mutation/synonym_replaced.csv")

head(one_hot$sample.ID)

head(sort_del)
