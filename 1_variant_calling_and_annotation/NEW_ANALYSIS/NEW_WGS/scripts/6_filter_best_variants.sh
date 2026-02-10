#!/bin/bash

#SBATCH -N 1
#SBATCH -n 1
#SBATCH -p medium
#SBATCH -J baseline_variant_set
#SBATCH -o /home/nikita.arya/GG20_WGS_NEW/log/baseline_variant_set.o%j
#SBATCH -e /home/nikita.arya/GG20_WGS_NEW/log/baseline_variant_set.e%j
#SBATCH --mail-user=nikita.a1@ahduni.edu.in
#SBATCH --mail-type=begin
#SBATCH --mail-type=end

module load anaconda3

source /home/nikita.arya/miniconda3/bin/activate gatk4

VARIANT_DIR="/home/nikita.arya/GG20_WGS_NEW/variant_calling"

gatk VariantFiltration \
   -V ${VARIANT_DIR}/raw_snps.vcf \
   -filter "QD < 15.0 || MQ < 60.00 || FS > 2.0 || SOR > 5.0 || MQRankSum < -3.09 || ReadPosRankSum < -1.7180" \
   --filter-name "best_quality_snps" \
   --missing-values-evaluate-as-failing \
   -O ${VARIANT_DIR}/baseline_snps.vcf

gatk SelectVariants \
   -V ${VARIANT_DIR}/baseline_snps.vcf \
   --exclude-filtered \
   -O ${VARIANT_DIR}/snps_pre_filtered.vcf

gatk VariantFiltration \
   -V ${VARIANT_DIR}/raw_indels.vcf \
   -filter "QD < 30.61 || MQ < 60.00 || FS > 0.01 || SOR > 0.37 || MQRankSum < -3.09 || ReadPosRankSum < -1.5440" \
   --filter-name "best_quality_indels" \
   --missing-values-evaluate-as-failing \
   -O ${VARIANT_DIR}/baseline_indels.vcf

gatk SelectVariants \
  -V ${VARIANT_DIR}/baseline_indels.vcf \
  --exclude-filtered \
  -O ${VARIANT_DIR}/indels_pre_filtered.vcf
