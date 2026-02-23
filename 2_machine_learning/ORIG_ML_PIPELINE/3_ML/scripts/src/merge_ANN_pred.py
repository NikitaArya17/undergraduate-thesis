import pandas as pd

df1 = pd.read_csv("3_ML/output_pred/pred_ANN_0_664.csv")
df2 = pd.read_csv("3_ML/output_pred/pred_ANN_664_1328.csv")
df3 = pd.read_csv("3_ML/output_pred/pred_ANN_1328_1991.csv")

df_list = [df1, df2, df3]

all_pred = pd.concat(df_list, axis = 0)

all_pred.to_csv("3_ML/output_pred/pred_ANN.csv", index = False)
