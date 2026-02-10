#!/bin/bash

#SBATCH -N 1
#SBATCH -n 8
#SBATCH -p medium
#SBATCH -J variant_filtration_tf
#SBATCH -o /home/nikita.arya/variant_filtration_tf.o%j
#SBATCH -e /home/nikita.arya/variant_filtration_tf.e%j
#SBATCH --mail-user=nikita.a1@ahduni.edu.in
#SBATCH --mail-type=begin
#SBATCH --mail-type=end

module load anaconda3

source /home/nikita.arya/miniconda3/bin/activate gatk4

# To filter SNPs
gatk VariantFiltration -R /home/nikita.arya/reference/Tifrunner/GCF_003086295_2_arahy_Tifrunner_gnm1_KYV3_genomic.fasta \
-V /home/nikita.arya/variants/Tifrunner/raw_snps_tf.vcf \
-filter-expression "(QD < 2.0) || (FS > 60.0) || (MQ < 40.0) || (SOR > 3.0)" \
-filter-name "basic_snp_filter" -O /home/nikita.arya/filtered_snps_tf.vcf

# To filter indels
gatk VariantFiltration -R /home/nikita.arya/reference/Tifrunner/GCF_003086295_2_arahy_Tifrunner_gnm1_KYV3_genomic.fasta \
-V /home/nikita.arya/variants/Tifrunner/raw_indels_tf.vcf \
-filter-expression "(QD < 2.0) || (FS > 60.0) || (MQ < 40.0) || (SOR > 3.0)" \
-filter-name "basic_indel_filter" -O /home/nikita.arya/filtered_indels_tf.vcf
