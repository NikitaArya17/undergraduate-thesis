#!/bin/bash

#SBATCH -N 1
#SBATCH -n 1
#SBATCH --cpus-per-task 8
#SBATCH -p long
#SBATCH -J bwa_alignment
#SBATCH -o /home/nikita.arya/GG20_WGS_NEW/log/bwa_alignment.o%j
#SBATCH -e /home/nikita.arya/GG20_WGS_NEW/log/bwa_alignment.e%j
#SBATCH --mail-user=nikita.a1@ahduni.edu.in
#SBATCH --mail-type=begin
#SBATCH --mail-type=end

module load anaconda3

source /home/nikita.arya/miniconda3/bin/activate gatk4

BASE_DIR="/home/nikita.arya/GG20_WGS_NEW"
READ1="${BASE_DIR}/trimmed_reads/51207407719-Stage-4-seeds-NCGM-4187_L001_R1_001_val_1.fq"
READ2="${BASE_DIR}/trimmed_reads/51207407719-Stage-4-seeds-NCGM-4187_L001_R2_001_val_2.fq"
REF_GENOME="${BASE_DIR}/ref_genome/ref_genome.fna"
ALIGN_DIR="${BASE_DIR}/alignment"

#INDEXING USING SAMTOOLS -
#samtools faidx ${REF_GENOME}

#INDEXING USING BWA-MEM2
bwa-mem2 index ${REF_GENOME}

bwa-mem2 mem -t 8 ${REF_GENOME} ${READ1} ${READ2} >> ${ALIGN_DIR}/aln.sam
