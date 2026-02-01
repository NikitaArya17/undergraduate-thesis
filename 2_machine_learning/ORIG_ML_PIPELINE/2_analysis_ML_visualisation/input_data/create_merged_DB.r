library(dplyr)

DB <- read.csv("1_create_DB/output_data/Database/DB_reproduction.csv")

str(DB)

merged_DB <- DB %>%
  arrange(culture_condition_ID)

str(merged_DB)

write.csv(merged_DB, "2_analysis_ML_visualisation/input_data/Merged_DB_Repro", row.names = FALSE)
