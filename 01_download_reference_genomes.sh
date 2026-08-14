#!/bin/bash
# =============================================================================
# SLURM JOB SUBMISSION: DOWNLOAD REFERENCE GENOMES
# Step 01 — run once, before submitting any downstream script. Downloads the
# two reference genomes this pipeline needs:
#   A) Chicken (Gallus gallus, GRCg7b) — used by 02_genome_assembly_array.sh
#      for RagTag pseudo-chromosome scaffolding and Liftoff annotation.
#   B) LEPC (Tympanuchus pallidicinctus, pur_lepc_1.0) — used by
#      05_combined_alignment_array.sh onward for the LEPC population-
#      genomics track (alignment, heterozygosity, PCA, ROH).
# Both downloads are independent and idempotent — skipped if already present.
# =============================================================================
#SBATCH --job-name=ref_genomes
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
# Chicken (Gallus gallus) RefSeq reference — GCF_016699485.2 / GRCg7b
CHICKEN_ASM_DIR="GCF_016699485.2_bGalGal1.mat.broiler.GRCg7b"
CHICKEN_FTP_BASE="https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/016/699/485/${CHICKEN_ASM_DIR}"
CHICKEN_REF_DIR="${CLUSTER_SCRATCH}/GROUSE_ASM/ref"
CHICKEN_FASTA="${CHICKEN_REF_DIR}/${CHICKEN_ASM_DIR}_genomic.fna"
CHICKEN_GFF="${CHICKEN_REF_DIR}/${CHICKEN_ASM_DIR}_genomic.gff"

# LEPC (Tympanuchus pallidicinctus) RefSeq reference — GCF_026119805.1 /
# pur_lepc_1.0. Output naming (.fa, not .fna; both compressed and
# decompressed kept) matches what 05/07/08 already expect — don't change
# without updating those scripts too.
LEPC_ASM_DIR="GCF_026119805.1_pur_lepc_1.0"
LEPC_FTP_BASE="https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/026/119/805/${LEPC_ASM_DIR}"
LEPC_REF_DIR="${CLUSTER_SCRATCH}/LEPC/ref"
LEPC_FASTA_GZ="${LEPC_REF_DIR}/${LEPC_ASM_DIR}_genomic.fa.gz"
LEPC_FASTA="${LEPC_REF_DIR}/${LEPC_ASM_DIR}_genomic.fa"

mkdir -p logs "$CHICKEN_REF_DIR" "$LEPC_REF_DIR"

echo ">>> 01_download_reference_genomes.sh"
echo ">>> Start time: $(date)"

# =============================================================================
# PART A: Chicken reference genome + GFF3 annotation
# =============================================================================
echo ">>> Part A: Chicken reference (${CHICKEN_ASM_DIR})"

if [[ -f "$CHICKEN_FASTA" ]]; then
    echo "    ${CHICKEN_FASTA} already present — skipping download"
else
    echo "    Downloading chicken reference genome"
    wget -q -O "${CHICKEN_FASTA}.gz" "${CHICKEN_FTP_BASE}/${CHICKEN_ASM_DIR}_genomic.fna.gz"
    gunzip -f "${CHICKEN_FASTA}.gz"
    echo "    Download complete: ${CHICKEN_FASTA}"
fi

if [[ -f "$CHICKEN_GFF" ]]; then
    echo "    ${CHICKEN_GFF} already present — skipping download"
else
    echo "    Downloading chicken annotation"
    wget -q -O "${CHICKEN_GFF}.gz" "${CHICKEN_FTP_BASE}/${CHICKEN_ASM_DIR}_genomic.gff.gz"
    gunzip -f "${CHICKEN_GFF}.gz"
    echo "    Download complete: ${CHICKEN_GFF}"
fi

# =============================================================================
# PART B: LEPC reference genome (no annotation needed downstream)
# =============================================================================
echo ">>> Part B: LEPC reference (${LEPC_ASM_DIR})"

if [[ -f "$LEPC_FASTA" ]]; then
    echo "    ${LEPC_FASTA} already present — skipping download"
else
    if [[ ! -f "$LEPC_FASTA_GZ" ]]; then
        echo "    Downloading LEPC reference genome"
        wget -q -O "$LEPC_FASTA_GZ" "${LEPC_FTP_BASE}/${LEPC_ASM_DIR}_genomic.fna.gz"
    fi
    echo "    Decompressing reference"
    gunzip -k -c "$LEPC_FASTA_GZ" > "$LEPC_FASTA"
    echo "    Download complete: ${LEPC_FASTA}"
fi

# =============================================================================
# NEXT STEP
# =============================================================================
echo ""
echo ">>> Step 01 complete."
echo ""
echo ">>> For the 3-species assembly pipeline, submit the assembly array with:"
echo "    N=\$(grep -vc '^#' assembly_manifest.tsv)"
echo "    sbatch --array=0-\$((N-2))%6 02_genome_assembly_array.sh"
echo ""
echo ">>> For the LEPC population-genomics pipeline, submit the combined"
echo ">>> alignment array's prep step with:"
echo "    sbatch --array=0-0 05_combined_alignment_array.sh --build-list-only"
echo ""
echo ">>> End time: $(date)"
