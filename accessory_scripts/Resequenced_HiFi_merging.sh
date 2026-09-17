#!/bin/bash
#SBATCH --job-name=grouse_hifi
#SBATCH -A fnrdewoody
#SBATCH -t 15:00:00
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --mem=40G
#SBATCH -p cpu
#SBATCH --mail-type=BEGIN,END,FAIL
#SBATCH --mail-user=blackan@purdue.edu

ml bioconda
ml anaconda
conda activate pbtk

#pbmerge -o F578_MERGED.bam F5478_m84221_260820_001659_s4.hifi_reads.bc2174.bam F5478_m84221_260901_214100_s2.hifi_reads.bc2174.bam

#pbmerge -o F5480_MERGED.bam F5480_pool_10_m84221_260821_081306_s3.hifi_reads.bc2175.bam F5480_pool_9_m84221_260821_060307_s2.hifi_reads.bc2175.bam

#pbmerge -o F5502_MERGED.bam F5502_m84221_260821_081306_s3.hifi_reads.bc2177.bam F5502_m84221_260901_214100_s2.hifi_reads.bc2177.bam

pbmerge -o F5595_MERGED.bam F5595_m84221_260819_220655_s3.hifi_reads.bc2186.bam F5595_m84221_260901_214100_s2.hifi_reads.bc2328.bam

mv F5595_m84221_* reseq/
