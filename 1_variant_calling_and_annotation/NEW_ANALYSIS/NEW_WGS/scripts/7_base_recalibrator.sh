#!/bin/bash

#SBATCH -N 1
#SBATCH -n 1
#SBATCH --mem=16G
#SBATCH -p long
#SBATCH -J base_recalibrator
#SBATCH -o /home/nikita.arya/GG20_WGS_NEW/log/base_recalibrator.o%j
#SBATCH -e /home/nikita.arya/GG20_WGS_NEW/log/base_recalibrator.e%j
#SBATCH --mail-user=nikita.a1@ahduni.edu.in
#SBATCH --mail-type=begin
#SBATCH --mail-type=end

module load anaconda3

source /home/nikita.arya/miniconda3/bin/activate gatk4

BASE_DIR="/home/nikita.arya/GG20_WGS_NEW"
REF_GENOME="${BASE_DIR}/ref_genome/ref_genome.fna"
ALIGN_DIR="${BASE_DIR}/alignment"
VARIANT_DIR="${BASE_DIR}/variant_calling"
BQSR_DIR="${BASE_DIR}/bqsr"
QUALI_DIR="${BASE_DIR}/qualimap"

gatk BaseRecalibrator -R ${REF_GENOME} \
-I ${ALIGN_DIR}/aln_dedup_sorted.bam \
--known-sites ${VARIANT_DIR}/snps_truth_set.vcf \
--known-sites ${VARIANT_DIR}/indels_truth_set.vcf \
-O ${BQSR_DIR}/recal.table

gatk ApplyBQSR -R ${REF_GENOME} \
-I ${ALIGN_DIR}/aln_dedup_sorted.bam \
-O ${BQSR_DIR}/recal_BSQR.bam -bqsr ${BQSR_DIR}/recal.table

gatk BaseRecalibrator -R ${REF_GENOME} \
-I ${BQSR_DIR}/recal_BSQR.bam \
--known-sites ${VARIANT_DIR}/snps_truth_set.vcf \
--known-sites ${VARIANT_DIR}/indels_truth_set.vcf \
-O ${BQSR_DIR}/postrecal.table

gatk AnalyzeCovariates -before ${BQSR_DIR}/recal.table -after ${BQSR_DIR}/postrecal.table \
-plots ${BQSR_DIR}/recal_plots.pdf

## NOTE: You should create a qualimap folder##

qualimap bamqc --java-mem-size=16G -bam ${BQSR_DIR}/recal_BSQR.bam \
-sd -oc ${QUALI_DIR}/genomecov.txt \
-ip -c -outdir ${QUALI_DIR} -outfile ${QUALI_DIR}/qualimap.pdf -outformat PDF
