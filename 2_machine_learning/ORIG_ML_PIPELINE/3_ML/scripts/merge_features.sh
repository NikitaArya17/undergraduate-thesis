#!/bin/bash

#SBATCH -n 1
#SBATCH --cpus-per-task 1
#SBATCH -p short
#SBATCH --mem=2G
#SBATCH -J merge_features
#SBATCH -o /home/nikita.arya/ensemble_pipeline/log/merge_features.o%j
#SBATCH -e /home/nikita.arya/ensemble_pipeline/log/merge_features.e%j
#SBATCH --mail-user=nikita.a1@ahduni.edu.in
#SBATCH --mail-type=begin
#SBATCH --mail-type=end


source /home/nikita.arya/miniconda3/bin/activate r_env

Rscript /home/nikita.arya/ensemble_pipeline/src/merge_features.r ANN #This can be run for all 3 model types by changing the end parameter to NB, SVM or ANN.
