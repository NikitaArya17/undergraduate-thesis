library(dplyr)

NB_data <- read.csv("2_analysis_ML_visualisation/output_data/NB_data_REPRODUCTION.csv", fileEncoding = "latin1")

dim(NB_data)

NB_data <- NB_data %>%
  select(!contains("Pert"), -Tem_high, -Strain_Evolved.MG1655)

colnames(NB_data) <- gsub("ÿ", "", colnames(NB_data))

NB_data[NB_data == "Absence"] <- 0
NB_data[NB_data == "Presence"] <- 1

write.csv(NB_data, "3_ML/output_data/MLDB_repro.csv", row.names = FALSE)
