#!/bin/bash

#SBATCH -n 1
#SBATCH --cpus-per-task 16
#SBATCH -p gpu
#SBATCH --mem=64G
#SBATCH -J ANN_feature_selection
#SBATCH -o /home/nikita.arya/ensemble_pipeline/YEAST_DATA/log/ANN_feature_selection.o%j
#SBATCH -e /home/nikita.arya/ensemble_pipeline/YEAST_DATA/log/ANN_feature_selection.e%j
#SBATCH --mail-user=nikita.a1@ahduni.edu.in
#SBATCH --mail-type=begin
#SBATCH --mail-type=end


source /home/nikita.arya/miniconda3/bin/activate python_ml

PARALLEL_LIMIT=$SLURM_CPUS_PER_TASK
CMD_FILE="/home/nikita.arya/ensemble_pipeline/YEAST_DATA/task_list_$SLURM_JOB_ID.txt"

echo "Generating task list..." > $CMD_FILE

for i in {0..7}
do

 for j in {47..49}
  do
    echo "python /home/nikita.arya/ensemble_pipeline/YEAST_DATA/src/Keras_FFN.py $i $j" >> $CMD_FILE
  done

done

TOTAL_JOBS=$(wc -l < $CMD_FILE)

echo "Generated $TOTAL_JOBS tasks. Starting parallel execution ($PARALLEL_LIMIT at a time)..."

cat $CMD_FILE | xargs -P $PARALLEL_LIMIT -I {} bash -c {}

rm $CMD_FILE

echo "All tasks completed"
