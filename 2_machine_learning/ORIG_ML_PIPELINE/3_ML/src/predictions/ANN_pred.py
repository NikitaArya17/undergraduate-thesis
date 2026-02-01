import pandas as pd
import numpy as np
import os
import keras
import keras.backend as K
from keras.models import Sequential,Model
from keras.layers import Dense, Dropout,BatchNormalization,Input
from keras.optimizers import RMSprop
from keras.regularizers import l2,l1
from keras.optimizers import Adam

import os
import collections


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


output_dir = "/home/nikita.arya/ensemble_pipeline/output_models/"


for i in range(1991):
    model_path = os.path.join(output_dir, f"ANN{i}.weights.h5")
    if os.path.exists(model_path):
        continue

    print(f"Training model {i} ...")

    Y_target = Y.iloc[:, i].values
    current_model = create_model(num_input = X.shape[1])
    current_model.fit(X, Y_target, epochs = 100, validation_split = 0.2, verbose = 0)
    current_model.save_weights(model_path)

    K.clear_session()
    del current_model


preds = []
test = X[0:1,:]
test[0:1,:] = -0.5
test[0:1,mask] = 0.5

for i in range(1991):
    model_path = f"/home/nikita.arya/ensemble_pipeline/output_models/ANN{i}.weights.h5"
    model.load_weights(model_path)

    pred  = model.predict(test)
    pred = np.squeeze(np.array(pred))
    preds.append(pred)
preds = np.array(preds)


res = pd.DataFrame()
res["genome site"] = df.columns[80:]
res["probability"] = preds


res.to_csv("/home/nikita.arya/ensemble_pipeline/output_predictions/pred_ANN.csv", index = False)
