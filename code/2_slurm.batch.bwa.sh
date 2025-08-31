#!/bin/bash

#SBATCH -N 1
#SBATCH -n 8
#SBATCH -p gpu
#SBATCH -J bwa_alignment
#SBATCH -o /home/nikita.arya/bwa_alignment.o%j
#SBATCH -e /home/nikita.arya/bwa_alignment.e%j
#SBATCH --mail-user=nikita.a1@ahduni.edu.in
#SBATCH --mail-type=begin
#SBATCH --mail-type=end

module load anaconda3

source /home/nikita.arya/miniconda3/bin/activate gatk4

#INDEXING USING SAMTOOLS - 
samtools faidx /home/nikita.arya/reference/GCF_003086295_2_arahy_Tifrunner_gnm1_KYV3_genomic.fasta

#INDEXING USING BWA-MEM2
bwa-mem2 index /home/nikita.arya/reference/GCF_003086295_2_arahy_Tifrunner_gnm1_KYV3_genomic.fasta

bwa-mem2 mem -t 8 /home/nikita.arya/reference/GCF_003086295_2_arahy_Tifrunner_gnm1_KYV3_genomic.fasta  \
/home/nikita.arya/raw_reads/GG20-leaf_R1.fastq.gz /home/nikita.arya/raw_reads/GG20-leaf_R2.fastq.gz >> /home/nikita.arya/alignment/aln.sam



