library(readxl)
library(dplyr)
library(stringr)
library(data.table)

phenos <- read_excel("YEAST_DATA/44320_2025_136_MOESM3_ESM.xlsx", sheet = 2)
genos_PAF <- read_excel("YEAST_DATA/44320_2025_136_MOESM11_ESM.xlsx", sheet = 2)
genos_PA <- fread("YEAST_DATA/genesMatrix_PresenceAbsence.tab", header = TRUE)


colnames(genos_PA)[1] <- colnames(genos_PAF)[1] <- "Standard_name"

phenos$Standard_name <- str_replace(phenos$Standard_name, "SACE_", "")
genos_PA$Standard_name <- str_replace(genos_PA$Standard_name, "SACE_", "")

phenos <- phenos %>%
  relocate(Standard_name)

phenos <- phenos[, -c(2:11)]
colnames(genos_PA)[1:20]
colnames(genos_PAF)[1:20]
colnames(phenos)[1:20]

df <- genos_PA %>%
  inner_join(genos_PAF, by = "Standard_name") %>%
  inner_join(phenos, by = "Standard_name")

write.csv(df, "YEAST_DATA/MLDB_yeast.csv", row.names = FALSE)
head(df)
dim(df)
class(df)

