#!/bin/bash

#SBATCH -N 1
#SBATCH -n 8
#SBATCH -p gpu
#SBATCH -J mark_duplicates
#SBATCH -o /home/nikita.arya/mark_duplicates.o%j
#SBATCH -e /home/nikita.arya/mark_duplicates.e%j
#SBATCH --mail-user=nikita.a1@ahduni.edu.in
#SBATCH --mail-type=begin
#SBATCH --mail-type=end

module load anaconda3

source /home/nikita.arya/miniconda3/bin/activate gatk4

#Be sure the file is sorted first; MarkDuplicates only works on sorted files
samtools sort -n /home/nikita.arya/alignment/aln.sam > /home/nikita.arya/sorted_dedup_alignment/aln.sorted.bam

gatk MarkDuplicates -O  /home/nikita.arya/sorted_dedup_alignment/aln_dedup.bam \
-M  /home/nikita.arya/metrics.txt \
-I  /home/nikita.arya/sorted_dedup_alignment/aln.sorted.bam --CREATE_INDEX true