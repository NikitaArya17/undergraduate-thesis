import pandas as pd
import openpyxl
from collections import defaultdict
from sklearn.preprocessing import MinMaxScaler

# 1) Scaling the input data to avoid biases in the ML predictions

pheno_index = pd.read_excel("YEAST_DATA/44320_2025_136_MOESM3_ESM.xlsx", sheet_name = 2)
all_phenos = pd.read_csv("YEAST_DATA/yeast_complete_phenotypes.csv", index_col = 0)

first_pheno = all_phenos.columns.get_loc('Sporulation_in_water')
first_pheno

gene_cols = list(all_phenos.columns[:first_pheno]) 

len(gene_cols)

unique_counts  = all_phenos[gene_cols].nunique()
valid_gene_cols = unique_counts[unique_counts > 1].index.tolist()

X = all_phenos[valid_gene_cols]
all_phenos[valid_gene_cols] = X.fillna(X.median())

scaler = MinMaxScaler()
all_phenos[valid_gene_cols] = scaler.fit_transform(all_phenos[valid_gene_cols])

#2) Converting the phenotypic data into categories suitable for classification

pheno_index.info()

categories = list(set(pheno_index['Defined clusters (Phenotype class)']))
categories

zipped_cols = zip(pheno_index['Defined clusters (Phenotype class)'], pheno_index['Phenotypes'])
paired_cols = list(zipped_cols)

trait_dict = defaultdict(list)

for category, trait in paired_cols:
    trait_dict[category].append(trait)

trait_dict = dict(trait_dict)

all_phenos.info()
all_phenos[categories] = None

for name in categories:
    valid_cols = [col for col in trait_dict.get(name, [])
                  if col in all_phenos.columns]
    
    if valid_cols:
        avg_score = all_phenos[valid_cols].mean(axis = 1)
        valid_mask = avg_score.notna()
        all_phenos.loc[valid_mask, name] = pd.qcut(avg_score.loc[valid_mask], 
                                                  q=2, labels=[0, 1]).astype(int)
    
    else:
        all_phenos[name] = None
        print(f'Warning: No columns found for category {name}')

all_phenos[categories].head()  

final_DB = all_phenos[valid_gene_cols + categories]
final_DB.shape

final_DB.to_csv("YEAST_DATA/MLDB_yeast.csv")
