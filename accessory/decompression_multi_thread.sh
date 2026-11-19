#!/bin/bash
#SBATCH --job-name=grouse_nexus
#SBATCH -A dewoody
#SBATCH -t 15:00:00
#SBATCH --nodes=1
#SBATCH --cpus-per-task=16
#SBATCH --mem=40G
#SBATCH -p cpu
#SBATCH --mail-type=BEGIN,END,FAIL
#SBATCH --mail-user=blackan@purdue.edu


THREADS=${SLURM_CPUS_PER_TASK}

tar -cvf ref_09_09_2026.tar.gz -I "pigz -p ${THREADS}" assembly/
