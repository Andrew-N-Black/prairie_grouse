#!/bin/bash
# =============================================================================
# SLURM ARRAY JOB: HI-C SUPPORT FOR EVERY RAGTAG JOIN (PTARMIGAN REFERENCE)
# Step 11 — requires
# 05_ragtag_liftoff_array.sh to have run with DO_HIC=true, which
# produces both the RagTag AGP and qc_ptarmigan/<sample>.<hap>.hic2final.cram.
#
# One array task per SAMPLE; each handles both haplotypes, matching steps 09
# and 10.
#
# THIS REUSES 11_join_support.py UNCHANGED — same scoring, same calibration,
# same thresholds. That is the point: the only thing that differs between the
# chicken run and this one is the reference that proposed the joins, so the two
# sets of scores are directly comparable. Keep 11_join_support.py in the
# submission directory.
#
# !! THE CRAM MUST BE THE PTARMIGAN ONE !!
#   The chicken-era qc/<sample>.<hap>.hic2final.cram is aligned to the
#   chicken-ORDERED assembly. The ptarmigan assembly holds the same sequence
#   under different names in a different order, so those coordinates do not
#   transfer and would yield confident nonsense. This script reads only from
#   qc_ptarmigan/, and step 05's DO_HIC stage is what fills it. If the per-
#   haplotype CRAM is missing the sample is skipped with a note rather than
#   scored against the wrong alignment.
#
# WHY THIS EXISTS
#   RagTag assembles each pseudo-chromosome by concatenating yahs scaffolds in
#   the order the REFERENCE implies, writing a 100 bp "align_genus" gap at every
#   join. Against chicken, 67 joins came out unsupported across the panel —
#   including a join at 19.1-19.5 Mb on chr_4 that is unsupported in 28
#   assemblies spanning all three species, which is where chicken's GGA4 fusion
#   sits, and which the literature says is a separate microchromosome in turkey
#   (MGA4 + MGA9).
#
#   This run asks whether those joins were telling us about the grouse or about
#   chicken. Reading the result:
#     - a chicken-unsupported join that is ABSENT here (ptarmigan never proposed
#       it) or SUPPORTED here  -> the chicken karyotype was the problem, and the
#       ptarmigan-ordered assembly is the better one to submit
#     - a join unsupported against BOTH references -> ours to explain: a real
#       assembly problem, or a rearrangement in Tympanuchus relative to both
#
#   Compare the two tables directly rather than reading this one alone; a join
#   count is only meaningful against its own denominator, which the controls
#   section of each run reports.
#
# HOW A JOIN IS SCORED (details in 11_join_support.py)
#   support = (observed - background) / (expected - background)
#     observed   contacts between flanking windows either side of the join
#     expected   what contiguous sequence gives at the same separations, from a
#                decay curve built only from WITHIN-scaffold bin pairs
#     background inter-chromosomal contact density
#   ~1.0 = as well supported as genuinely contiguous sequence. ~0.0 = no more
#   linked than two different chromosomes.
#
#   The script measures its own calibration: it scores pseudo-joins at interior
#   positions of large intact scaffolds, which are contiguous by construction
#   and must land near 1.0. Verdict thresholds are set as fractions of that
#   control median, so they adapt to each library rather than relying on a
#   number tuned elsewhere. If the controls do not centre sensibly the script
#   withholds verdicts instead of reporting numbers it cannot stand behind.
#
# USAGE
#   N=$(grep -v '^#' assembly_manifest.tsv | tail -n +2 | grep -c .)
#   sbatch --array=0-$((N-1))%8 11_join_support.sh
#
# COLLATING AFTERWARDS
#   cd $CLUSTER_SCRATCH/GROUSE/grouse_asm/qc_ptarmigan/join_support
#   head -1 $(ls *.joins.tsv | head -1) > ALL.joins.tsv
#   tail -q -n +2 *.joins.tsv >> ALL.joins.tsv
#   awk -F'\t' 'NR>1 && $12=="unsupported"' ALL.joins.tsv | sort -k11,11g | head -40
# =============================================================================
#SBATCH --job-name=grouse_joinsupport_ptarmigan
#SBATCH --output=logs/%x_%A_%a.out
#SBATCH --error=logs/%x_%A_%a.err
#SBATCH -A dewoody
#SBATCH -t 08:00:00
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
# 24G: the binned contact matrix dominates — at the default 100 kb bins a
# ~1.05 Gb assembly gives ~10,500 bins, so a dense int32 matrix is ~440 MB.
# The rest is samtools decode buffers. Haplotypes are processed one after the
# other, so the peak is per-haplotype rather than additive.
#SBATCH --mem=24G
#SBATCH -p cpu
#SBATCH --mail-type=BEGIN,END,FAIL
#SBATCH --mail-user=blackan@purdue.edu

set -euo pipefail

module unload anaconda 2>/dev/null || true
ml biocontainers
ml samtools/1.22.1

# unset LD_PRELOAD: RCAC's XALT library is injected this way and fails on some
# nodes (GLIBC mismatch), which can kill subshells under `set -e`.
unset LD_PRELOAD

# =============================================================================
# USER-DEFINED VARIABLES
# =============================================================================
MANIFEST="${SLURM_SUBMIT_DIR}/assembly_manifest.tsv"

PROJECT_DIR="${CLUSTER_SCRATCH}/GROUSE/grouse_asm"

# Must match REF_TAG in 05_ragtag_liftoff_array.sh. Every input path
# is derived from it, which is what keeps this run from reading the chicken-era
# CRAM (see the warning in the header).
REF_TAG="ptarmigan"

FINAL_DIR="${PROJECT_DIR}/final_${REF_TAG}"
RAGTAG_DIR="${PROJECT_DIR}/ragtag_${REF_TAG}"
QC_DIR="${PROJECT_DIR}/qc_${REF_TAG}"
OUT_DIR="${QC_DIR}/join_support"

SCORE_SCRIPT="${SLURM_SUBMIT_DIR}/11_join_support.py"

# ---- Scoring parameters -----------------------------------------------------
# Contact-matrix bin size. 100 kb balances resolution at a join against having
# enough read pairs per bin. Raise it if the script reports the decay curve as
# too sparse to sit above background.
BIN_SIZE=100000

# Flank size either side of a join, in bins. 10 bins x 100 kb = 1 Mb per side.
WINDOW_BINS=10

# Hi-C reads below this MAPQ are dropped, matching the assembly pipeline.
MIN_MAPQ=20

# Only score joins where BOTH flanking scaffolds are at least this long. Short
# scaffolds carry too few contacts for the statistic to mean anything, and
# joining a 70 kb fragment is low-stakes either way.
MIN_COMPONENT=1000000

THREADS=$SLURM_CPUS_PER_TASK

mkdir -p logs "$OUT_DIR"

# =============================================================================
# PRE-FLIGHT
# =============================================================================
if [[ ! -f "$SCORE_SCRIPT" ]]; then
    echo "ERROR: companion script not found: ${SCORE_SCRIPT}"
    echo "       11_join_support.py must sit next to this file."
    exit 1
fi
for BIN in samtools python3; do
    if ! command -v "$BIN" >/dev/null 2>&1; then
        echo "ERROR: '${BIN}' not on PATH."
        exit 1
    fi
done

# numpy is the only non-stdlib dependency; reuse the env step 09 built if the
# system python3 lacks it.
PY="$(command -v python3)"
if ! "$PY" -c "import numpy" >/dev/null 2>&1; then
    CAND="${PROJECT_DIR}/conda_envs/buscoplot/bin/python3"
    if [[ -x "$CAND" ]] && "$CAND" -c "import numpy" >/dev/null 2>&1; then
        PY="$CAND"
    else
        echo "ERROR: no python3 with numpy available."
        echo "       Run 09_busco_plots.sh once to create the plotting env, or"
        echo "       load a python module that provides numpy."
        exit 1
    fi
fi
echo ">>> python   : ${PY}  [$("$PY" -c 'import numpy;print("numpy",numpy.__version__)')]"
echo ">>> samtools : $(command -v samtools)"

if [[ -z "${SLURM_ARRAY_TASK_ID:-}" ]]; then
    echo "ERROR: SLURM_ARRAY_TASK_ID is not set."
    echo "Submit with: sbatch --array=0-N 11_join_support.sh"
    exit 1
fi
if [[ ! -f "$MANIFEST" ]]; then
    echo "ERROR: Manifest not found: ${MANIFEST}"
    exit 1
fi

mapfile -t ROWS < <(grep -v '^#' "$MANIFEST" | tail -n +2 | grep -v '^[[:space:]]*$')
LINE="${ROWS[$SLURM_ARRAY_TASK_ID]:-}"
if [[ -z "$LINE" ]]; then
    echo "ERROR: No manifest row at index ${SLURM_ARRAY_TASK_ID} (${#ROWS[@]} samples)"
    exit 1
fi
IFS=$'\t' read -r SAMPLE SPECIES _REST <<< "$LINE"
if [[ -z "${SAMPLE:-}" || -z "${SPECIES:-}" ]]; then
    echo "ERROR: Malformed manifest row: ${LINE}"
    exit 1
fi

echo ">>> Array task ${SLURM_ARRAY_TASK_ID} -> ${SAMPLE} (${SPECIES})"
echo ">>> Started: $(date)"

# =============================================================================
# PER-HAPLOTYPE
# =============================================================================
for HAP in hap1 hap2; do
    PREFIX="${SPECIES}_${SAMPLE}_${HAP}"
    FINAL_FASTA="${FINAL_DIR}/${PREFIX}.pseudo_chr.fasta"
    FAI="${FINAL_FASTA}.fai"
    AGP="${RAGTAG_DIR}/${SAMPLE}.${HAP}/ragtag.scaffold.agp"
    RENAME="${RAGTAG_DIR}/${SAMPLE}.${HAP}/${SAMPLE}.${HAP}.rename_map.tsv"
    CRAM="${QC_DIR}/${SAMPLE}.${HAP}.hic2final.cram"
    OUT="${OUT_DIR}/${PREFIX}"

    echo ""
    echo "============================================================"
    echo ">>> ${PREFIX}"
    echo "============================================================"

    MISSING=false
    for F in "$FINAL_FASTA" "$FAI" "$AGP" "$CRAM"; do
        if [[ ! -s "$F" ]]; then
            echo "  ! missing: ${F}"
            MISSING=true
        fi
    done
    if [[ "$MISSING" == true ]]; then
        echo "  skipping ${PREFIX} — step 05 outputs incomplete"
        echo "    If only the .hic2final.cram is missing, step 05 ran with"
        echo "    DO_HIC=false. Set it true and rerun step 05 for this sample;"
        echo "    the chicken-era CRAM under qc/ is NOT a substitute."
        continue
    fi

    if [[ -s "${OUT}.joins.tsv" && "${OUT}.joins.tsv" -nt "$CRAM" ]]; then
        echo "  already scored — skipping"
        continue
    fi

    # Stream the Hi-C alignment once. -F 0x904 drops unmapped, secondary and
    # supplementary records up front so they never reach python.
    samtools view -@ "$THREADS" -F 0x904 -T "$FINAL_FASTA" "$CRAM" \
        | "$PY" "$SCORE_SCRIPT" \
            --agp "$AGP" \
            --fai "$FAI" \
            --rename-map "$RENAME" \
            --out-prefix "$OUT" \
            --assembly "$PREFIX" \
            --bin-size "$BIN_SIZE" \
            --window-bins "$WINDOW_BINS" \
            --min-mapq "$MIN_MAPQ" \
            --min-component "$MIN_COMPONENT"
done

echo ""
echo "============================================================"
echo ">>> ${SAMPLE} complete: $(date)"
echo "  ${OUT_DIR}/${SPECIES}_${SAMPLE}_hap{1,2}.joins.tsv"
echo "  ${OUT_DIR}/${SPECIES}_${SAMPLE}_hap{1,2}.controls.tsv"
echo "============================================================"
