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
zcat 20260901_1622_3A_PBM69309_f3a9da84_Sample_F5468/fastq_pass_F5468/*fastq.gz > F5468.fastq ; gzip F5468.fastq
seqkit stats 20260901_1622_3A_PBM69309_f3a9da84_Sample_F5468/fastq_pass_F5468/*fastq.gz -o  multi_fastqs_F5468.txt
seqkit stats F5468.fastq.gz -o single_fastq_F5468.txt
seqkit rmdup -n F5468.fastq.gz -o  F5468_dedup.fastq.gz -d duplicated_ids_F5468.txt

zcat 20260901_1519_3B_PBM67548_2762ccfd_Sample_F5545/fastq_pass_F5545/*fastq.gz > F5545.fastq ; gzip F5545.fastq
seqkit stats 20260901_1519_3B_PBM67548_2762ccfd_Sample_F5545/fastq_pass_F5545/*fastq.gz  -o  multi_fastqs_F5545.txt
seqkit stats F5545.fastq.gz -o single_fastqs_F5545.txt
seqkit rmdup -n F5545.fastq.gz -o  F5545_dedup.fastq.gz -d duplicated_ids_F5545.txt

zcat 20260901_1622_3C_PBM66906_d8a57dea_Sample_F5597/fastq_pass_F5597/*fastq.gz > F5597.fastq ; gzip F5597.fastq
seqkit stats 20260901_1622_3C_PBM66906_d8a57dea_Sample_F5597/fastq_pass_F5597/*fastq.gz  -o multi_fastqs_F5597.txt
seqkit stats F5597.fastq.gz -o single_fastqs_F5597.txt
seqkit rmdup -n F5597.fastq.gz -o  F5597_dedup.fastq.gz -d duplicated_ids_F5597.txt
