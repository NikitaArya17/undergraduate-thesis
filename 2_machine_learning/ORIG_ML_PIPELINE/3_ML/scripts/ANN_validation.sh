#!/bin/bash

#SBATCH -n 1
#SBATCH --cpus-per-task 8
#SBATCH -p long
#SBATCH --mem=32G
#SBATCH -J ANN_validation
#SBATCH -o /home/nikita.arya/ensemble_pipeline/log/ANN_validation.o%j
#SBATCH -e /home/nikita.arya/ensemble_pipeline/log/ANN_validation.e%j
#SBATCH --mail-user=nikita.a1@ahduni.edu.in
#SBATCH --mail-type=begin
#SBATCH --mail-type=end


source /home/nikita.arya/miniconda3/bin/activate python_ml


python /home/nikita.arya/ensemble_pipeline/src/FFN_val.py
