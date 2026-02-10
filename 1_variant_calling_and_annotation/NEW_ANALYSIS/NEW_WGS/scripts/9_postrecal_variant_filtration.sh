#!/bin/bash

#SBATCH -N 1
#SBATCH -n 1
#SBATCH -p long
#SBATCH -J variant_refiltration
#SBATCH -o /home/nikita.arya/GG20_WGS_NEW/log/variant_refiltration.o%j
#SBATCH -e /home/nikita.arya/GG20_WGS_NEW/log/variant_refiltration.e%j
#SBATCH --mail-user=nikita.a1@ahduni.edu.in
#SBATCH --mail-type=begin
#SBATCH --mail-type=end

module load anaconda3

source /home/nikita.arya/miniconda3/bin/activate gatk4

BASE_DIR="/home/nikita.arya/GG20_WGS_NEW"
REF_GENOME="${BASE_DIR}/ref_genome/ref_genome.fna"
VARIANT_DIR="${BASE_DIR}/variant_calling"

# To filter SNPs
gatk VariantFiltration -R ${REF_GENOME} \
-filter-expression "QD < 31.75 || MQ < 60.00 || FS > 0.01 || SOR > 0.37 || MQRankSum < -5.33 || ReadPosRankSum < -1.7180" \
-filter-name "final_snp_filter" \
-V ${VARIANT_DIR}/raw_recal_snps.vcf \
-O ${VARIANT_DIR}/filtered_recal_snps.vcf

#To filter indels
#gatk VariantFiltration -R ${REF_GENOME} \
#-filter-expression "QD < 30.61 || MQ < 60.00 || FS > 0.01 || SOR > 0.37 || MQRankSum < -3.09 || ReadPosRankSum < -1.5440" \
#-filter-name "final_indel_filter" \
#-V ${VARIANT_DIR}/raw_recal_indels.vcf \
#-O ${VARIANT_DIR}/filtered_recal_indels.vcf

# Selecting the unfiltered variants only
#gatk SelectVariants -V ${VARIANT_DIR}/filtered_recal_indels.vcf -O ${VARIANT_DIR}/filtered_only_indels.vcf \
#--select "vc.isNotFiltered()"

gatk SelectVariants -V ${VARIANT_DIR}/filtered_recal_snps.vcf -O ${VARIANT_DIR}/filtered_only_snps.vcf \
--select "vc.isNotFiltered()"
