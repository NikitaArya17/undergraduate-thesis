#!/bin/bash

#SBATCH -N 1
#SBATCH -n 1
#SBATCH --cpus-per-task 1
#SBATCH --mem=8G
#SBATCH -p long
#SBATCH -J filter_variant_tranches
#SBATCH -o /home/nikita.arya/GG20_WGS_NEW/log/filter_variant_tranches.o%j
#SBATCH -e /home/nikita.arya/GG20_WGS_NEW/log/filter_variant_tranches.e%j
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

gatk FilterVariantTranches \
   -V ${VARIANT_DIR}/recal_annotated_variants.vcf \
   --resource ${VARIANT_DIR}/snps_truth_set.vcf \
   --resource ${VARIANT_DIR}/indels_truth_set.vcf \
   --info-key CNN_2D \
   --snp-tranche 95.0 --snp-tranche 99.0 \
   --indel-tranche 90.0 --indel-tranche 95.0 \
   --invalidate-previous-filters \
   -O ${VARIANT_DIR}/final_filtered_variants.vcf

gatk SelectVariants -R ${REF_GENOME} \
-V ${VARIANT_DIR}/final_filtered_variants.vcf -select-type SNP \
-O ${VARIANT_DIR}/final_filtered_snps.vcf

gatk SelectVariants -R ${REF_GENOME} \
-V ${VARIANT_DIR}/final_filtered_variants.vcf -select-type INDEL \
-O ${VARIANT_DIR}/final_filtered_indels.vcf
