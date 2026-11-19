#!/bin/bash
#SBATCH --job-name=grouse_ONT
#SBATCH -A dewoody
#SBATCH -t 15:00:00
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --mem=80G
#SBATCH -p cpu
#SBATCH --mail-type=BEGIN,END,FAIL
#SBATCH --mail-user=blackan@purdue.edu

ml biocontainers seqkit

seqkit stats F5597.fastq.gz -o single_fastqs_F5597.txt

seqkit seq -m 100000 F5468.fastq.gz > F5468_ul_100kb.fastq.gz
seqkit seq -m 100000 F5545.fastq.gz > F5545_ul_100kb.fastq.gz
seqkit seq -m 100000 F5597.fastq.gz > F5597_ul_100kb.fastq.gz

seqkit stats -a F5468_ul_100kb.fastq.gz -o F5468_ul_100kb_stats.txt
seqkit stats -a F5545_ul_100kb.fastq.gz -o F5545_ul_100kb_stats.txt
seqkit stats -a  F5597_ul_100kb.fastq.gz -o F5597_ul_100kb_stats.txt
