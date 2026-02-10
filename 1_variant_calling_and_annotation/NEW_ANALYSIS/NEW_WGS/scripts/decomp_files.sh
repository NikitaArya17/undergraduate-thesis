#!/bin/bash

#SBATCH -N 1
#SBATCH -n 1
#SBATCH -p short
#SBATCH -J decomp_files
#SBATCH -o /home/nikita.arya/GG20_WGS_NEW/log/decomp_files.o%j
#SBATCH -e /home/nikita.arya/GG20_WGS_NEW/log/decomp_files.e%j

BASE_DIR="/home/nikita.arya/GG20_WGS_NEW"

gunzip ${BASE_DIR}/raw_reads/*
