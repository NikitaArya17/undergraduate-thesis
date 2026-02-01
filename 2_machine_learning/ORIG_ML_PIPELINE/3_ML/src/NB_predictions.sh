#!/bin/bash

#SBATCH -n 1
#SBATCH --cpus-per-task 1
#SBATCH -p medium
#SBATCH --mem=4G
#SBATCH -J NB_predictions
#SBATCH -o /home/nikita.arya/ensemble_pipeline/log/NB_predictions.o%j
#SBATCH -e /home/nikita.arya/ensemble_pipeline/log/NB_predictions.e%j
#SBATCH --mail-user=nikita.a1@ahduni.edu.in
#SBATCH --mail-type=begin
#SBATCH --mail-type=end


source /home/nikita.arya/miniconda3/bin/activate r_env

Rscript /home/nikita.arya/ensemble_pipeline/src/NB_pred.R
