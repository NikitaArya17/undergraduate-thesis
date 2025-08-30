#!/bin/bash

#SBATCH -N 1
#SBATCH -n 8
#SBATCH -p short
#SBATCH -J fast_qc
#SBATCH -o /home/nikita.arya/fast_qc.o%j
#SBATCH -e /home/nikita.arya/fast_qc.e%j
#SBATCH --mail-user=nikita.a1@ahduni.edu.in
#SBATCH --mail-type=begin
#SBATCH --mail-type=end

module load anaconda3

source /home/nikita.arya/miniconda3/bin/activate gatk4

fastqc -o /home/nikita.arya/fastqc_output/ -t 8 \
/home/nikita.arya/raw_reads/GG20-leaf_R1.fastq.gz /home/nikita.arya/raw_reads/GG20-leaf_R2.fastq.gz
