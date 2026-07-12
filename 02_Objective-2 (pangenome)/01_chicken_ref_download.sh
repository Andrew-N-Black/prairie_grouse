#!/bin/bash
# =============================================================================
# SLURM JOB SUBMISSION: DOWNLOAD CHICKEN REFERENCE (Gallus gallus, GRCg7b)
# Step 01 — run once, before submitting 02_genome_assembly_array.sh.
# Downloads the reference FASTA + GFF3 annotation used by that array job's
# RagTag pseudo-chromosome scaffolding and Liftoff annotation steps.
# =============================================================================
#SBATCH --job-name=chicken_ref
#SBATCH --output=logs/%x_%j.out
#SBATCH --error=logs/%x_%j.err
#SBATCH -A fnrdewoody
#SBATCH -t 0-04:00:00
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=8G
#SBATCH -p cpu
#SBATCH --mail-type=BEGIN,END,FAIL
#SBATCH --mail-user=${USER}@purdue.edu

# =============================================================================
# ENVIRONMENT SETUP
# =============================================================================
set -euo pipefail

# =============================================================================
# USER-DEFINED VARIABLES
# =============================================================================
PROJECT="GROUSE_ASM"
PROJECT_DIR="${CLUSTER_SCRATCH}/${PROJECT}"
REF_DIR="${PROJECT_DIR}/ref"

# Chicken (Gallus gallus) RefSeq reference — GCF_016699485.2 / GRCg7b
CHICKEN_ACC="GCF_016699485.2"
CHICKEN_ASM_NAME="bGalGal1.mat.broiler.GRCg7b"
CHICKEN_ASM_DIR="${CHICKEN_ACC}_${CHICKEN_ASM_NAME}"
CHICKEN_FTP_BASE="https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/016/699/485/${CHICKEN_ASM_DIR}"
CHICKEN_FASTA="${REF_DIR}/${CHICKEN_ASM_DIR}_genomic.fna"
CHICKEN_GFF="${REF_DIR}/${CHICKEN_ASM_DIR}_genomic.gff"

mkdir -p logs "$REF_DIR"

echo ">>> 02_download_chicken_ref.sh"
echo ">>> Reference : ${CHICKEN_ACC} (${CHICKEN_ASM_NAME})"
echo ">>> Start time: $(date)"

# =============================================================================
# STEP 1: Download reference genome (skipped if already present)
# =============================================================================
if [[ -f "$CHICKEN_FASTA" ]]; then
    echo ">>> ${CHICKEN_FASTA} already present — skipping download"
else
    echo ">>> Downloading chicken reference genome"
    wget -q -O "${CHICKEN_FASTA}.gz" "${CHICKEN_FTP_BASE}/${CHICKEN_ASM_DIR}_genomic.fna.gz"
    gunzip -f "${CHICKEN_FASTA}.gz"
    echo "    Download complete: ${CHICKEN_FASTA}"
fi

# =============================================================================
# STEP 2: Download GFF3 annotation (skipped if already present)
# =============================================================================
if [[ -f "$CHICKEN_GFF" ]]; then
    echo ">>> ${CHICKEN_GFF} already present — skipping download"
else
    echo ">>> Downloading chicken annotation"
    wget -q -O "${CHICKEN_GFF}.gz" "${CHICKEN_FTP_BASE}/${CHICKEN_ASM_DIR}_genomic.gff.gz"
    gunzip -f "${CHICKEN_GFF}.gz"
    echo "    Download complete: ${CHICKEN_GFF}"
fi

# =============================================================================
# NEXT STEP
# =============================================================================
echo ""
echo ">>> Step 02 complete. Submit the assembly array with:"
echo ""
echo "    N=\$(grep -vc '^#' assembly_manifest.tsv)"
echo "    sbatch --array=0-\$((N-2))%6 03_genome_assembly_array.sh"
echo ""
echo ">>> End time: $(date)"
