#!/bin/bash
# =============================================================================
# accessory/te_categories.sh -- TE composition by category (CR1, ERVs, ...) per haplotype, for the autosomes,
# the Z and the W, from the RepeatMasker output of 12_pangenome_build.sh (STAGE=te).
#   sbatch accessory/te_categories.sh        # from the repo root
# Output: ${COMPARE_DIR}/te_categories.tsv
# =============================================================================
#SBATCH -A dewoody
#SBATCH -p cpu
#SBATCH -q standby
#SBATCH -t 04:00:00
#SBATCH -c 1
#SBATCH --mem=8G
#SBATCH -J te_cat
#SBATCH -o logs/%x_%j.out
set -euo pipefail
REPO_DIR="${SLURM_SUBMIT_DIR:-$PWD}"
PROJECT_DIR="${CLUSTER_SCRATCH}/GROUSE/grouse_asm"
FINAL_DIR="${FINAL_DIR:-${PROJECT_DIR}/final_ptarmigan}"
TE_DIR="${TE_DIR:-${PROJECT_DIR}/te_annotation}"
COMPARE_DIR="${COMPARE_DIR:-${PROJECT_DIR}/comparison}"
MANIFEST="${MANIFEST:-${REPO_DIR}/assembly_manifest.tsv}"
python3 "${REPO_DIR}/accessory/te_categories.py" "$MANIFEST" "$FINAL_DIR" "$TE_DIR" "${COMPARE_DIR}/te_categories.tsv"
echo ">>> ${COMPARE_DIR}/te_categories.tsv"
