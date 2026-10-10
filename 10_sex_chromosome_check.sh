#!/bin/bash
# =============================================================================
# SLURM JOB: IS THE HAPLOTYPE-LENGTH ASYMMETRY A MISSING Z CHROMOSOME?
# Step 10 — requires 07_busco_array.sh. Check 3 additionally requires
# 08_ptarmigan_reference_busco.sh; without it checks 1 and 2 still run and
# check 3 is skipped with a note.
#
# Six assemblies (F5595 hap2, F5596 hap2, F5598 hap1, F5599 hap1, F5600 hap1
# and F5503 hap1) are both ~100 Mb shorter than their partner haplotype AND
# ~4 percentage points lower in BUSCO completeness. This script tests, rather
# than infers, whether the shortfall is loss of Z-linked sequence.
#
#   1. chr_Z / chr_W / chr_MT length per assembly vs the panel median
#   2. how many Complete BUSCOs each assembly places on its own chr_Z
#   3. where each assembly's MISSING orthologs sit on the ptarmigan reference,
#      and what fraction of the excess missing is Z-linked -> verdict
#
# Checks 1 and 2 are independent of the reference run and of each other, so a
# consistent answer across all three is strong evidence; a split answer tells
# you which assumption broke.
#
# This is seconds of work (parsing .fai and full_table.tsv, pure standard
# library), so it runs just as happily on a login node:
#
#     bash 10_sex_chromosome_check.sh
#
# or as a job:
#
#     sbatch 10_sex_chromosome_check.sh
#
# OUTPUT, under ${QC_DIR}/sex_chromosome_check:
#   chromosome_content.tsv            per assembly: total, autosome, Z, W, MT
#   busco_by_chromosome.tsv           per assembly: Complete BUSCOs on Z and W
#   missing_busco_localisation.tsv    per assembly: missing orthologs by
#                                     reference chromosome (check 3 only)
# =============================================================================
#SBATCH --job-name=grouse_sexcheck_ptarmigan
#SBATCH --output=logs/%x_%j.out
#SBATCH --error=logs/%x_%j.err
#SBATCH -A fnrdewoody
#SBATCH -t 00:30:00
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
# 8G: the largest structure in memory is 47 BUSCO full tables held as dicts of
# ~8,338 entries each. Single-threaded and I/O-bound.
#SBATCH --mem=8G
#SBATCH -p cpu
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=blackan@purdue.edu

set -euo pipefail

# unset LD_PRELOAD: RCAC's XALT library is injected this way and fails on some
# nodes (GLIBC mismatch), which can kill subshells under `set -e`.
unset LD_PRELOAD

# =============================================================================
# USER-DEFINED VARIABLES
# =============================================================================
PROJECT_DIR="${CLUSTER_SCRATCH}/GROUSE/grouse_asm"

# Must match REF_TAG in 05_ragtag_liftoff_array.sh.
REF_TAG="ptarmigan"

FINAL_DIR="${PROJECT_DIR}/final_${REF_TAG}"
QC_DIR="${PROJECT_DIR}/qc_${REF_TAG}"
BUSCO_DIR="${QC_DIR}/busco"
OUT_DIR="${QC_DIR}/sex_chromosome_check"

# Directory name of the ptarmigan reference BUSCO run, as written by
# 08_ptarmigan_reference_busco.sh.
REF_NAME="LAGMUT_bLagMut1_ref"

# An assembly is called Z-deficient when its chr_Z is below this fraction of
# the panel median chr_Z. 0.5 is deliberately loose: a haplotype that has lost
# the Z retains only the pseudoautosomal and gametologous fraction, which in
# F5503 hap1 was 9.4 Mb against roughly 76 Mb, i.e. about 0.12x.
Z_DEFICIT_FRAC=0.5

# Works from the submission directory when run with sbatch, and from the
# current directory when run as a plain bash script on a login node.
SCRIPT_DIR="${SLURM_SUBMIT_DIR:-$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)}"
CHECK_SCRIPT="${SCRIPT_DIR}/10_sex_chromosome_check.py"

mkdir -p logs "$OUT_DIR"

# =============================================================================
# PRE-FLIGHT
# =============================================================================
if [[ ! -f "$CHECK_SCRIPT" ]]; then
    echo "ERROR: companion script not found: ${CHECK_SCRIPT}"
    echo "       10_sex_chromosome_check.py must sit next to this file."
    exit 1
fi
if [[ ! -d "$BUSCO_DIR" ]]; then
    echo "ERROR: BUSCO output directory not found: ${BUSCO_DIR}"
    echo "       Run 07_busco_array.sh first."
    exit 1
fi

# The analysis is pure standard library, so any python3 will do. Prefer the
# env built in step 09 when it exists, purely because its version is known.
PLOT_PYTHON="${PROJECT_DIR}/conda_envs/buscoplot/bin/python3"
if [[ -x "$PLOT_PYTHON" ]]; then
    PY="$PLOT_PYTHON"
elif command -v python3 >/dev/null 2>&1; then
    PY="$(command -v python3)"
else
    echo "ERROR: no python3 found."
    exit 1
fi
echo ">>> python: ${PY}  [$("$PY" --version 2>&1)]"

# =============================================================================
# RUN
# =============================================================================
echo ">>> Started: $(date)"
echo ""

"$PY" "$CHECK_SCRIPT" \
    --final-dir "$FINAL_DIR" \
    --busco-dir "$BUSCO_DIR" \
    --out-dir   "$OUT_DIR" \
    --ref-name  "$REF_NAME" \
    --z-deficit-frac "$Z_DEFICIT_FRAC"

echo ""
echo "============================================================"
echo ">>> Complete: $(date)"
echo "  ${OUT_DIR}/chromosome_content.tsv"
echo "  ${OUT_DIR}/busco_by_chromosome.tsv"
echo "  ${OUT_DIR}/missing_busco_localisation.tsv   (check 3 only)"
echo "============================================================"
