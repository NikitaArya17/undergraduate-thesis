#!/bin/bash

#SBATCH -N 1
#SBATCH -n 8
#SBATCH -p gpu
#SBATCH -J samtools_stats
#SBATCH -o /home/nikita.arya/samtools_stats.o%j
#SBATCH -e /home/nikita.arya/samtools_stats.e%j
#SBATCH --mail-user=nikita.a1@ahduni.edu.in
#SBATCH --mail-type=begin
#SBATCH --mail-type=end

module load anaconda3

source /home/nikita.arya/miniconda3/bin/activate gatk4

# Collecting depth, read alignment quality and other metrics
#Sort the file by genomic coordinates first
samtools sort /home/nikita.arya/sorted_dedup_alignment/aln_dedup.bam > /home/nikita.arya/sorted_dedup_alignment/aln_dedup_sorted.bam

samtools flagstat /home/nikita.arya/sorted_dedup_alignment/aln_dedup_sorted.bam > /home/nikita.arya/GG20_stats.txt

samtools stats /home/nikita.arya/sorted_dedup_alignment/aln_dedup_sorted.bam > /home/nikita.arya/GG20_reads.stats

gatk CollectAlignmentSummaryMetrics -R /home/nikita.arya/reference/GCF_003086295_2_arahy_Tifrunner_gnm1_KYV3_genomic.fasta \
-I /home/nikita.arya/sorted_dedup_alignment/aln_dedup_sorted.bam \
-O /home/nikita.arya/alignmentmetrics.txt

gatk CollectInsertSizeMetrics -I /home/nikita.arya/sorted_dedup_alignment/aln_dedup_sorted.bam \
-O /home/nikita.arya/insertsize.txt \
-H /home/nikita.arya/insertsize.pdf