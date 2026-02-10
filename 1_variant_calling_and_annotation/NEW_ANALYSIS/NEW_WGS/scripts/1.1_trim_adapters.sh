#!/bin/bash
#SBATCH --job-name=trim_adapters       # Job name
#SBATCH --partition=medium             # Specify the partition/queue
#SBATCH --nodes=1                    # Run on a single node
#SBATCH --cpus-per-task=4            # Number of cores
#SBATCH --output=/home/nikita.arya/GG20_WGS_NEW/log/trim_adapters.o%j
#SBATCH --error=/home/nikita.arya/GG20_WGS_NEW/log/trim_adapters.e%j
#SBATCH --mail-user=nikita.a1@ahduni.edu.in
#SBATCH --mail-type=BEGIN
#SBATCH --mail-type=END

# --- SETUP --- #
# Load the main anaconda module for the HPC environment
module load anaconda3

# Define the main working directory and output directories
BASE_DIR="/home/nikita.arya/GG20_WGS_NEW"
OUTPUT_DIR="${BASE_DIR}/trimmed_reads"

source /home/nikita.arya/miniconda3/bin/activate trim_adapters

# Define file names for the current sample
READ1="${BASE_DIR}/raw_reads/51207407719-Stage-4-seeds-NCGM-4187_L001_R1_001.fastq"
READ2="${BASE_DIR}/raw_reads/51207407719-Stage-4-seeds-NCGM-4187_L001_R2_001.fastq"

trim_galore --paired -q 20 --length 20 -j 4 --fastqc --output_dir ${OUTPUT_DIR} ${READ1} ${READ2}

echo "--- Adapter trimming complete for sample: ${CURRENT_SAMPLE} ---"
