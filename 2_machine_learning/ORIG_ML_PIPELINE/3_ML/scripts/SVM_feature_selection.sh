#!/bin/bash

#SBATCH -n 1
#SBATCH --cpus-per-task 10
#SBATCH -p long
#SBATCH --mem=32G
#SBATCH -J SVM_feature_selection
#SBATCH -o /home/nikita.arya/ensemble_pipeline/log/SVM_feature_selection.o%j
#SBATCH -e /home/nikita.arya/ensemble_pipeline/log/SVM_feature_selection.e%j
#SBATCH --mail-user=nikita.a1@ahduni.edu.in
#SBATCH --mail-type=begin
#SBATCH --mail-type=end


source /home/nikita.arya/miniconda3/bin/activate r_env

PARALLEL_LIMIT=$SLURM_CPUS_PER_TASK
CMD_FILE="task_list_$SLURM_JOB_ID.txt"

echo "Generating task list..." > $CMD_FILE

for i in {1..1990}
do
    echo "Rscript /home/nikita.arya/ensemble_pipeline/src/feature_selection.R $i 2" >> $CMD_FILE
done

TOTAL_JOBS=$(wc -l < $CMD_FILE)

echo "Generated $TOTAL_JOBS tasks. Starting parallel execution ($PARALLEL_LIMIT at a time)..."

cat $CMD_FILE | xargs -P $PARALLEL_LIMIT -I {} bash -c {}

rm $CMD_FILE

echo "All tasks completed"
