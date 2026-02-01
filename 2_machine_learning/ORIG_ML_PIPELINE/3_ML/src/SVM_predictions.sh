#!/bin/bash

#SBATCH -n 1
#SBATCH --cpus-per-task 2
#SBATCH -p long
#SBATCH --mem=8G
#SBATCH -J SVM_predictions
#SBATCH -o /home/nikita.arya/ensemble_pipeline/log/SVM_predictions.o%j
#SBATCH -e /home/nikita.arya/ensemble_pipeline/log/SVM_predictions.e%j
#SBATCH --mail-user=nikita.a1@ahduni.edu.in
#SBATCH --mail-type=begin
#SBATCH --mail-type=end


source /home/nikita.arya/miniconda3/bin/activate r_env

Rscript /home/nikita.arya/ensemble_pipeline/src/SVM_pred.R
