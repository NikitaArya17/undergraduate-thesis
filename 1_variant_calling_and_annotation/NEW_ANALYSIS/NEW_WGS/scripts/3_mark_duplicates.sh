#!/bin/bash

#SBATCH -N 1
#SBATCH -n 1
#SBATCH --cpus-per-task 6
#SBATCH -p long
#SBATCH -J mark_duplicates
#SBATCH -o /home/nikita.arya/GG20_WGS_NEW/log/mark_duplicates.o%j
#SBATCH -e /home/nikita.arya/GG20_WGS_NEW/log/mark_duplicates.e%j
#SBATCH --mail-user=nikita.a1@ahduni.edu.in
#SBATCH --mail-type=begin
#SBATCH --mail-type=end

module load anaconda3

source /home/nikita.arya/miniconda3/bin/activate gatk4

BASE_DIR="/home/nikita.arya/GG20_WGS_NEW"
ALIGN_DIR="${BASE_DIR}/alignment"
METRICS_DIR="${BASE_DIR}/stats_metrics"

#Sort the input file
samtools sort -@ 6 -o ${ALIGN_DIR}/aln.sorted.bam ${ALIGN_DIR}/aln.sam

#Add Read Groups to the file and sort the output
gatk AddOrReplaceReadGroups \
    -I ${ALIGN_DIR}/aln.sorted.bam \
    -O ${ALIGN_DIR}/aln.rg.sorted.bam \
    --RGID ARNDAVIIL5.1 \
    --RGLB roche_prep_01 \
    --RGPL ILLUMINA \
    --RGPU ARNDAVIIL5.1.AGATACGG+CCGTCCTT \
    --RGSM GG20_reads \
    --SORT_ORDER coordinate \
    --CREATE_INDEX true

gatk MarkDuplicates -O  ${ALIGN_DIR}/aln_dedup.bam \
-M  ${METRICS_DIR}/dedup_metrics.txt \
-I  ${ALIGN_DIR}/aln.rg.sorted.bam --CREATE_INDEX true
