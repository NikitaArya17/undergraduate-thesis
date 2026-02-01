import pandas as pd
import numpy as np
import os
import tensorflow as tf
from tensorflow.keras.models import Sequential
from tensorflow.keras.layers import Dense
from tensorflow.keras.optimizers import Adam
from sklearn.metrics import roc_auc_score, accuracy_score


base_dir = "/home/nikita.arya/ensemble_pipeline/YEAST_DATA"
input_file = os.path.join(base_dir, "input_data/MLDB_yeast.csv")
feature_dir = os.path.join(base_dir, "output_data/feature_selection/ANN")


pred_output_file = os.path.join(base_dir, "output_predictions/ANN/ANN_validation.csv")
metrics_output_file = os.path.join(base_dir, "output_predictions/ANN/ANN_metrics_summary.csv")


os.makedirs(os.path.dirname(pred_output_file), exist_ok=True)


print("Loading Yeast Data...")
df = pd.read_csv(input_file, index_col=0)
pheno_names = df.columns[-8:].tolist()
X_raw = df.iloc[:, :-8].values


pred_df = pd.DataFrame(index=df.index)
metrics_list = [] # We will store AUC/accuracy here


for p in pheno_names:
    pred_df[f"Actual_{p}"] = df[p]


def build_model(input_dim):
    model = Sequential()
    model.add(Dense(32, input_dim=input_dim, activation='relu'))
    model.add(Dense(16, activation='relu'))
    model.add(Dense(1, activation='sigmoid'))
    model.compile(loss='binary_crossentropy', optimizer=Adam(learning_rate=0.001), metrics=['accuracy'])
    return model


for pheno in pheno_names:
    print(f"\nProcessing {pheno}...")


    y_target = df[pheno].values
    mask = ~np.isnan(y_target)

    for setting in [47, 48, 49]:
        feat_path = os.path.join(feature_dir, f"pheno_class_{pheno_names.index(pheno)}_setting_{setting}")


        if not os.path.exists(feat_path):
            continue

        try:
            with open(feat_path, 'r') as f:

                line = f.readline().replace(',', ' ')
                indices = [int(float(x)) for x in line.split() if x.replace('.', '', 1).isdigit()]
        except:
            print(f"  Error reading {feat_path}")
            continue

        if not indices: continue


        X_sub = X_raw[:, indices]


        X_train = X_sub[mask]
        y_train = y_target[mask]


        model = build_model(len(indices))

        model.fit(X_train, y_train, epochs=150, batch_size=16, verbose=0)



        probs = model.predict(X_sub, verbose=0).flatten()


        col_name = f"Pred_{pheno}_Set_{setting}"
        pred_df[col_name] = probs


        val_probs = probs[mask]

        try:
            auc = roc_auc_score(y_train, val_probs)
            # Threshold at 0.5 for accuracy
            acc = accuracy_score(y_train, (val_probs > 0.5).astype(int))
        except ValueError:
            auc = 0.5
            acc = 0.0

        print(f"  -> Set {setting} ({len(indices)} genes) | AUC: {auc:.3f}")


        metrics_list.append({
            'Phenotype': pheno,
            'Setting': setting,
            'Num_Genes': len(indices),
            'AUC': auc,
            'Accuracy': acc
        })

print("\nSaving raw predictions...")
pred_df.to_csv(pred_output_file)

print("Saving metrics summary...")
metrics_df = pd.DataFrame(metrics_list)
metrics_df.to_csv(metrics_output_file, index=False)

print(f"Done! Files saved to:\n1. {pred_output_file}\n2. {metrics_output_file}")
