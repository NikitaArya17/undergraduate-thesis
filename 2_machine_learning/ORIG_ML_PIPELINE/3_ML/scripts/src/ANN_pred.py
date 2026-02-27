import pandas as pd
import numpy as np
import gc
import keras
from sklearn.model_selection import train_test_split
import keras.backend as K
from keras.models import Sequential,Model
from keras.layers import Dense, Dropout,BatchNormalization,Input
from keras.optimizers import RMSprop
from keras.regularizers import l2,l1
from keras.optimizers import Adam

import sys
import os
import collections

os.environ["CUDA_VISIBLE_DEVICES"] = "-1"

start_idx = 0
end_idx = 1991

if len(sys.argv) > 2:
    start_idx = int(sys.argv[1])
    end_idx = int(sys.argv[2])
    print(f"Processing genes {start_idx} to {end_idx}")

df = pd.read_csv("/home/nikita.arya/ensemble_pipeline/output_data/MLDB_repro.csv")
setting = pd.read_csv("/home/nikita.arya/ensemble_pipeline/output_data/CultureCondition.csv", dtype={'Temperature':object})

setting = ["_".join(map(str, pair)) for pair in zip(setting.columns,setting.iloc[0,:])]
mask = df.columns.isin(setting)[:81]

first_gene_index = df.columns.get_loc("ybcS")

X, Y = np.split(df, [first_gene_index], axis=1)
X = X.values
X = X-0.5

Model_setting = collections.namedtuple('Model_setting','num_layers num_node alpha drop_rate act_method lr regularization patience')


setting_ = [2,100, 0.5, 0.2, 'tanh', 0.01, 'l2', 3]
setting = Model_setting(*setting_)
setting = setting._asdict()


num_output_ = 1
def create_model(num_input = 79,num_output = num_output_):
    X_input = Input(shape=(num_input,))

    X = Dense(64)(X_input)
    X = Dropout(0.2)(X)
    X = Dense(32)(X)
    Ys= []
    for i in range(num_output):
        Ys.append(Dense(1, activation = 'sigmoid')(X))
    model = Model(inputs=[X_input],outputs = Ys)
    model.compile(loss=['binary_crossentropy']*num_output,loss_weights=[1.]*num_output,optimizer=Adam(learning_rate=setting['lr']), metrics=['accuracy'])
    return model


output_dir = "/home/nikita.arya/ensemble_pipeline/output_models/ANN"

num_samples = X.shape[0]
indices = np.arange(num_samples)
train_idx, val_idx = train_test_split(indices, test_size=0.2, random_state=2021)

for i in range(start_idx, end_idx):
    model_path = os.path.join(output_dir, f"ANN{i}.weights.h5")
    if os.path.exists(model_path):
        continue

    print(f"Training model {i} ...")

    Y_target = Y.iloc[:, i].values

    X_train, X_val = X[train_idx], X[val_idx]
    y_train, y_val = Y_target[train_idx], Y_target[val_idx]

    current_model = create_model(num_input = X.shape[1])
    current_model.fit(X_train, y_train,
                      epochs=100,
                      validation_data=(X_val, y_val),
                      verbose=0)
    current_model.save_weights(model_path)

    K.clear_session()
    del current_model
    gc.collect()


preds = []
test = X[0:1,:].copy()
test[0:1,:] = -0.5
test[0:1,mask] = 0.5

for i in range(start_idx, end_idx):
    model = create_model(num_input = X.shape[1])
    model_path = f"/home/nikita.arya/ensemble_pipeline/output_models/ANN/ANN{i}.weights.h5"
    model.load_weights(model_path)

    pred  = model.predict(test, verbose = 0)
    pred = np.squeeze(np.array(pred))
    preds.append(pred)
    K.clear_session()
    del model
    gc.collect()

preds = np.array(preds)


res = pd.DataFrame()
all_genes = df.columns[first_gene_index:]
res["genome_site"] = all_genes[start_idx:end_idx]
res["probability"] = preds

output_file = f"/home/nikita.arya/ensemble_pipeline/output_predictions/ANN/pred_ANN_{start_idx}_{end_idx}.csv"
res.to_csv(output_file, index = False)
print(f"Saved predictions to {output_file}")
