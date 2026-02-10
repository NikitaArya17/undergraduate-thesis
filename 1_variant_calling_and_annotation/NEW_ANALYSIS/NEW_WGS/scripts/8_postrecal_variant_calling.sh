#!/bin/bash

#SBATCH -N 1
#SBATCH -n 1
#SBATCH --cpus-per-task 4
#SBATCH --mem=16G
#SBATCH -p long
#SBATCH -J haplotype_caller
#SBATCH -o /home/nikita.arya/GG20_WGS_NEW/log/haplotype_caller.o%j
#SBATCH -e /home/nikita.arya/GG20_WGS_NEW/log/haplotype_caller.e%j
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

#Sort the input file
samtools sort -@ 4 -o ${BQSR_DIR}/recal_BQSR_sorted.bam ${BQSR_DIR}/recal_BSQR.bam

#Index the input file first
samtools index ${BQSR_DIR}/recal_BQSR_sorted.bam

gatk HaplotypeCaller -ploidy 4 -R ${REF_GENOME} \
-I ${BQSR_DIR}/recal_BQSR_sorted.bam \
-O ${VARIANT_DIR}/raw_recal_variants.vcf

gatk SelectVariants -R ${REF_GENOME} \
-V ${VARIANT_DIR}/raw_recal_variants.vcf -select-type SNP \
-O ${VARIANT_DIR}/raw_recal_snps.vcf

gatk SelectVariants -R ${REF_GENOME} \
-V ${VARIANT_DIR}/raw_recal_variants.vcf -select-type INDEL \
-O ${VARIANT_DIR}/raw_recal_indels.vcf
