#!/bin/bash

#SBATCH -N 1
#SBATCH -n 8
#SBATCH -p medium
#SBATCH -J base_recalibrator_tf
#SBATCH -o /home/nikita.arya/base_recalibrator_tf.o%j
#SBATCH -e /home/nikita.arya/base_recalibrator_tf.e%j
#SBATCH --mail-user=nikita.a1@ahduni.edu.in
#SBATCH --mail-type=begin
#SBATCH --mail-type=end

module load anaconda3

source /home/nikita.arya/miniconda3/bin/activate gatk4

gatk BaseRecalibrator -R /home/nikita.arya/reference/Tifrunner/GCF_003086295_2_arahy_Tifrunner_gnm1_KYV3_genomic.fasta \
-I /home/nikita.arya/sorted_dedup_alignment/Tifrunner/aln_dedup_sorted.bam \
--known-sites /home/nikita.arya/variants/Tifrunner/filtered_snps_tf.vcf \
--known-sites /home/nikita.arya/variants/Tifrunner/filtered_indels_tf.vcf \
-O /home/nikita.arya/recal_tf.table

gatk ApplyBQSR -R /home/nikita.arya/reference/Tifrunner/GCF_003086295_2_arahy_Tifrunner_gnm1_KYV3_genomic.fasta \
-I /home/nikita.arya/sorted_dedup_alignment/Tifrunner/aln_dedup_sorted.bam \
-O /home/nikita.arya/recal_BSQR_tf.bam -bqsr /home/nikita.arya/recal_tf.table

gatk BaseRecalibrator -R /home/nikita.arya/reference/Tifrunner/GCF_003086295_2_arahy_Tifrunner_gnm1_KYV3_genomic.fasta \
-I /home/nikita.arya/recal_BSQR_tf.bam \
--known-sites /home/nikita.arya/variants/Tifrunner/filtered_snps_tf.vcf \
--known-sites /home/nikita.arya/variants/Tifrunner/filtered_indels_tf.vcf \
-O /home/nikita.arya/postrecal_tf.table

gatk AnalyzeCovariates -before /home/nikita.arya/recal_tf.table -after /home/nikita.arya/postrecal_tf.table \
-plots /home/nikita.arya/recal_plots_tf.pdf

## NOTE: You should create a qualimap folder##

qualimap bamqc --java-mem-size=200G -bam /home/nikita.arya/recal_BSQR_tf.bam \
-sd -oc /home/nikita.arya/genomecov_tf.txt \
-ip -c -outdir /home/nikita.arya/qualimap/Tifrunner -outfile /home/nikita.arya/qualimap/Tifrunner/qualimap_tf.pdf -outformat PDF
