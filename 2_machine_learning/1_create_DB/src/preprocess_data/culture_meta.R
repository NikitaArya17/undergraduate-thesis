source("1_create_DB/src/preprocess_data/setpath.R") 
source("1_create_DB/src/helper/culture_condition/culture_one_hot_encode.R")

input_path = file.path(data_path,"culture_condition.csv")
output_path = file.path(out_data_path,"1_culture_matrix.csv")

culture_one_hot_encode(input_path,output_path) # corrected function name
