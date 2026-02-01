import pandas as pd
import numpy as np
import gc
import keras
import keras.backend as K
from keras.models import Sequential,Model
from keras.layers import Dense, Dropout,BatchNormalization,Input
from keras.optimizers import RMSprop
from keras.regularizers import l2,l1
from keras.optimizers import Adam
from sklearn.feature_selection import SelectKBest, f_classif
import os
import collections

base_dir = "/home/nikita.arya/ensemble_pipeline/YEAST_DATA"
output_dir = f"{base_dir}/output_models/ANN"

df = pd.read_csv(f"{base_dir}/input_data/MLDB_yeast.csv", index_col = 0)

first_pheno_index = df.shape[1] - 8
gene_features = df.columns[:first_pheno_index]
pheno_targets = df.columns[first_pheno_index:]


Model_setting = collections.namedtuple('Model_setting','num_layers num_node alpha drop_rate act_method lr regularization patience')

setting_ = [2,100, 0.5, 0.2, 'tanh', 0.01, 'l2', 3]
setting = Model_setting(*setting_)
setting = setting._asdict()


def create_model(num_input):
    model = Sequential()
    model.add(Dense(64, input_shape=(num_input,), activation=setting['act_method']))
    model.add(Dropout(setting['drop_rate']))
    model.add(Dense(32, activation=setting['act_method']))
    model.add(Dense(1, activation='sigmoid'))

    model.compile(
        optimizer=Adam(learning_rate=setting['lr']),
        loss='binary_crossentropy',
        metrics=['accuracy']
    )
    return model


selected_features_map = {}

for i, target_name in enumerate(pheno_targets):
    print(f"\n--- Processing Phenotype {i}: {target_name} ---")

    current_df = df[list(gene_features) + [target_name]].dropna(subset=[target_name])

    X_all = current_df[gene_features].values
    y_all = current_df[target_name].values

    selector = SelectKBest(score_func=f_classif, k=100)
    X_selected = selector.fit_transform(X_all, y_all)

    selected_indices = selector.get_support(indices=True)
    selected_features_map[target_name] = selected_indices

    model_path = os.path.join(output_dir, f"ANN_{target_name}.weights.h5")

    model = create_model(num_input=100)
    model.fit(
        X_selected, y_all,
        epochs=100,
        validation_split=0.2,
        verbose=0,
        batch_size=32
    )
    model.save_weights(model_path)

    print(f"Model saved for {target_name} using {X_selected.shape[1]} features.")

    K.clear_session()
    del model
    gc.collect()


print("\n--- Starting Prediction ---")
results_df = pd.DataFrame(index=df.index)

for target_name in pheno_targets:
    model = create_model(num_input=100)
    model_path = os.path.join(output_dir, f"ANN_{target_name}.weights.h5")
    model.load_weights(model_path)

    indices = selected_features_map[target_name]
    X_test = df[gene_features].iloc[:, indices].values


    probs = model.predict(X_test, verbose=0)
    results_df[f"Prob_{target_name}"] = probs.flatten()

    results_df[f"Actual_{target_name}"] = df[target_name]

    K.clear_session()
    del model
    gc.collect()

final_output = os.path.join(base_dir, "output_predictions", "ANN", "ANN_pred.csv")
results_df.to_csv(final_output)
print(f"\nSuccess! Final predictions for all strains saved to: {final_output}")
