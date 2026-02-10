#!/bin/bash

#SBATCH -N 1
#SBATCH -n 1
#SBATCH -p medium
#SBATCH -J variant_selection
#SBATCH --array=0-1
#SBATCH -o /home/nikita.arya/GG20_WGS_NEW/log/variant_selection_%A_%a.out
#SBATCH -e /home/nikita.arya/GG20_WGS_NEW/log/variant_selection_%A_%a.err
#SBATCH --mail-user=nikita.a1@ahduni.edu.in
#SBATCH --mail-type=begin
#SBATCH --mail-type=end

module load anaconda3

source /home/nikita.arya/miniconda3/bin/activate gatk4

BASE_DIR="/home/nikita.arya/GG20_WGS_NEW/variant_calling"
TYPES=("snps" "indels")
MUT_TYPE=${TYPES[$SLURM_ARRAY_TASK_ID]}

gatk VariantsToTable \
-V ${BASE_DIR}/raw_${MUT_TYPE}_tf.vcf \
-F CHROM -F POS -F QD -F MQ -F FS -F SOR -F MQRankSum -F ReadPosRankSum \
-O ${BASE_DIR}/raw_${MUT_TYPE}.table
