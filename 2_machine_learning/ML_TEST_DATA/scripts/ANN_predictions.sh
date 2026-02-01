#!/bin/bash

#SBATCH -n 1
#SBATCH --cpus-per-task 10
#SBATCH -p gpu
#SBATCH --mem=32G
#SBATCH -J ANN_predictions
#SBATCH -o /home/nikita.arya/ensemble_pipeline/YEAST_DATA/log/ANN_predictions.o%j
#SBATCH -e /home/nikita.arya/ensemble_pipeline/YEAST_DATA/log/ANN_predictions.e%j
#SBATCH --mail-user=nikita.a1@ahduni.edu.in
#SBATCH --mail-type=begin
#SBATCH --mail-type=end


source /home/nikita.arya/miniconda3/bin/activate python_ml

python /home/nikita.arya/ensemble_pipeline/YEAST_DATA/src/ANN_pred.py
