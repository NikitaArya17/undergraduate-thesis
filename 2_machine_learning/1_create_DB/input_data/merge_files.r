library(dplyr)
library(tidyr)
library(stringr)

cul_cond <- read.csv("1_create_DB/input_data/culture_condition.csv")
mut_uni <- read.csv("1_create_DB/input_data/mutation_before_unification.csv")

head(cul_cond)
colnames(cul_cond)

head(mut_uni)
colnames(mut_uni)

genes <- data.frame(mut_uni$Sample.ID, mut_uni$PMID, mut_uni$Mutation.type, mut_uni$Gene)

head(genes)

genes_sep_rows <- genes %>%
  separate_longer_delim(cols = "mut_uni.Gene", delim = ";")

head(genes_sep_rows)
dim(genes_sep_rows)

unique(genes_sep_rows$mut_uni.Mutation.type)

# for uniformity in mutation type names
genes_sep_rows$mut_uni.Mutation.type <- str_replace_all(genes_sep_rows$mut_uni.Mutation.type,
".*NP.*", "SNP")

genes_sep_rows$mut_uni.Mutation.type <- str_replace_all(genes_sep_rows$mut_uni.Mutation.type,
"SNP\n", "SNP")

genes_sep_rows$mut_uni.Mutation.type <- str_replace_all(genes_sep_rows$mut_uni.Mutation.type,
".*el.*", "Deletion")

genes_sep_rows$mut_uni.Mutation.type <- str_replace_all(genes_sep_rows$mut_uni.Mutation.type,
".*mp.*", "Amplification")

genes_sep_rows$mut_uni.Mutation.type <- str_replace_all(genes_sep_rows$mut_uni.Mutation.type,
".*ser.*", "Insertion")

testing <- cul_cond %>%
  inner_join(genes_sep_rows, by = c( "ID" = "mut_uni.Sample.ID", "PMID" = "mut_uni.PMID"))

## Try transposing half and rejoining

str(testing)
dim(testing)

colnames(testing)[32] <- "Mutation.type"
colnames(testing)[33] <- "Gene"

str(testing)
