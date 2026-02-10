#!/bin/bash

#SBATCH -N 1
#SBATCH -n 1
#SBATCH --cpus-per-task 2
#SBATCH -p medium
#SBATCH -J fast_qc
#SBATCH -o /home/nikita.arya/GG20_WGS_NEW/log/fast_qc.o%j
#SBATCH -e /home/nikita.arya/GG20_WGS_NEW/log/fast_qc.e%j
#SBATCH --mail-user=nikita.a1@ahduni.edu.in
#SBATCH --mail-type=begin
#SBATCH --mail-type=end

module load anaconda3

source /home/nikita.arya/miniconda3/bin/activate gatk4

fastqc -o /home/nikita.arya/GG20_WGS_NEW/output_data/fastqc_output/ -t 2 \
/home/nikita.arya/GG20_WGS_NEW/raw_reads/51207407719-Stage-4-seeds-NCGM-4187_L001_R1_001.fastq.gz \
/home/nikita.arya/GG20_WGS_NEW/raw_reads/51207407719-Stage-4-seeds-NCGM-4187_L001_R2_001.fastq.gz
