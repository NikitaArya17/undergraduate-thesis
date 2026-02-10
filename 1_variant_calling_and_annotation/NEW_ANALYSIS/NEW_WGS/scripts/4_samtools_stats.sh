#!/bin/bash

#SBATCH -N 1
#SBATCH -n 1
#SBATCH -p medium
#SBATCH -J samtools_stats
#SBATCH -o /home/nikita.arya/GG20_WGS_NEW/log/samtools_stats.o%j
#SBATCH -e /home/nikita.arya/GG20_WGS_NEW/log/samtools_stats.e%j
#SBATCH --mail-user=nikita.a1@ahduni.edu.in
#SBATCH --mail-type=begin
#SBATCH --mail-type=end

module load anaconda3

source /home/nikita.arya/miniconda3/bin/activate gatk4

BASE_DIR="/home/nikita.arya/GG20_WGS_NEW"
REF_GENOME="${BASE_DIR}/ref_genome/ref_genome.fna"
ALIGN_DIR="${BASE_DIR}/alignment"
METRICS_DIR="${BASE_DIR}/stats_metrics"

# Collecting depth, read alignment quality and other metrics

samtools flagstat ${ALIGN_DIR}/aln_dedup.bam > ${METRICS_DIR}/GG20_stats.txt

samtools stats ${ALIGN_DIR}/aln_dedup.bam > ${METRICS_DIR}/GG20_reads.stats

gatk CollectAlignmentSummaryMetrics -R ${REF_GENOME} \
-I ${ALIGN_DIR}/aln_dedup.bam \
-O ${METRICS_DIR}/alignmentmetrics.txt

gatk CollectInsertSizeMetrics -I ${ALIGN_DIR}/aln_dedup.bam \
-O ${METRICS_DIR}/insertsize.txt \
-H ${METRICS_DIR}/insertsize.pdf
