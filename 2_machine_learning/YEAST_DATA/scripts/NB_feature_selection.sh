#!/bin/bash

#SBATCH -n 1
#SBATCH --cpus-per-task 4
#SBATCH -p long
#SBATCH --mem=8G
#SBATCH -J NB_feature_selection
#SBATCH -o /home/nikita.arya/ensemble_pipeline/YEAST_DATA/log/NB_feature_selection.o%j
#SBATCH -e /home/nikita.arya/ensemble_pipeline/YEAST_DATA/log/NB_feature_selection.e%j
#SBATCH --mail-user=nikita.a1@ahduni.edu.in
#SBATCH --mail-type=begin
#SBATCH --mail-type=end

source /home/nikita.arya/miniconda3/bin/activate r_env

CMD_FILE="nb_task_list_$SLURM_JOB_ID.txt" > $CMD_FILE

for p in {0..7}
do
    for s in 47 48 49
    do
        echo "Rscript /home/nikita.arya/ensemble_pipeline/YEAST_DATA/src/NB_feature_selection.R $p $s" >> $CMD_FILE
    done
done

cat $CMD_FILE | xargs -P 4 -I {} bash -c "{}"

rm $CMD_FILE
echo "All Naive Bayes tasks completed."
