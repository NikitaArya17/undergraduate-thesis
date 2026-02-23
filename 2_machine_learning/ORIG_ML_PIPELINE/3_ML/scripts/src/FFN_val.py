import pandas as pd
import numpy as np
import os
import tensorflow as tf
from tensorflow.keras.models import Model
from tensorflow.keras.layers import Dense, Dropout, Input
from tensorflow.keras.optimizers import Adam
from tensorflow.keras import backend as K
from sklearn.model_selection import train_test_split
from sklearn.metrics import roc_auc_score, accuracy_score, confusion_matrix
import gc

base_dir = "/home/nikita.arya/ensemble_pipeline/"
input_file = os.path.join(base_dir, "output_data/MLDB_repro.csv")
model_dir = os.path.join(base_dir, "output_models/ANN")

output_dir = os.path.join(base_dir, "output_data/ANN_validation")
os.makedirs(output_dir, exist_ok=True)

pred_output_file = os.path.join(output_dir, "ANN_validation.csv")
metrics_output_file = os.path.join(output_dir, "ANN_metrics.csv")

print("Loading data...")
df = pd.read_csv(input_file)

if "ybcS" in df.columns:
    first_gene_index = df.columns.get_loc("ybcS")
else:
    first_gene_index = 80

print(f"Split index: {first_gene_index}")
print(f"Features: {first_gene_index} (Input Shape)")

X_all = df.iloc[:, :first_gene_index].values
X_all = X_all - 0.5

gene_names = df.columns[first_gene_index:]
Y_all = df.iloc[:, first_gene_index:].values

def create_model(num_input):
    X_input = Input(shape=(num_input,))
    X = Dense(64)(X_input)
    X = Dropout(0.2)(X)
    X = Dense(32)(X)

    Ys = []
    Ys.append(Dense(1, activation='sigmoid')(X))

    model = Model(inputs=[X_input], outputs=Ys)
    model.compile(loss='binary_crossentropy', optimizer='adam', metrics=['accuracy'])
    return model

metrics_list = []
val_preds_df = pd.DataFrame(np.nan, index=df.index, columns=gene_names)

indices = np.arange(len(df))
_, val_idx = train_test_split(indices, test_size=0.2, random_state=42)
X_val = X_all[val_idx]

print(f"Starting validation for {len(gene_names)} genes...")

for i, gene_name in enumerate(gene_names):
    weight_file = os.path.join(model_dir, f"ANN{i}.weights.h5")

    if not os.path.exists(weight_file):
        continue

    try:
        model = create_model(num_input=X_all.shape[1])
        model.load_weights(weight_file)


        preds = model.predict(X_val, verbose=0)

        if isinstance(preds, list):
            probs = preds[0].flatten()
        else:
            probs = preds.flatten()


        val_preds_df.iloc[val_idx, i] = probs


        y_true = Y_all[val_idx, i]

        try:
            auc = roc_auc_score(y_true, probs)
        except ValueError:
            auc = 0.5

        y_pred_class = (probs > 0.5).astype(int)
        acc = accuracy_score(y_true, y_pred_class)

        tn, fp, fn, tp = confusion_matrix(y_true, y_pred_class, labels=[0, 1]).ravel()
        sensitivity = tp / (tp + fn) if (tp + fn) > 0 else 0
        specificity = tn / (tn + fp) if (tn + fp) > 0 else 0

        metrics_list.append({
            'Gene_ID': i,
            'Gene_Name': gene_name,
            'AUC': auc,
            'Accuracy': acc,
            'Sensitivity': sensitivity,
            'Specificity': specificity
        })

        K.clear_session()
        del model

    except Exception as e:
        print(f"Error processing {gene_name} (Index {i}): {e}")
        continue

    if i % 100 == 0:
        print(f"Processed {i}/{len(gene_names)}")
        gc.collect()

print("\nSaving files...")
metrics_df = pd.DataFrame(metrics_list)
metrics_df.to_csv(metrics_output_file, index=False)
val_preds_df.to_csv(pred_output_file)

print("Done.")
