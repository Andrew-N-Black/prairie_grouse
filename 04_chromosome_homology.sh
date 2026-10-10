#!/bin/bash
# =============================================================================
# SLURM JOB: WHAT DOES EACH PTARMIGAN CHROMOSOME CORRESPOND TO IN CHICKEN?
# Step 04 — requires 03_ptarmigan_reference.sh for the ptarmigan reference and
# 01_download_chicken_reference.sh for the chicken one.
#
# This is the step that lets a ptarmigan chromosome number be related to the
# chicken karyotype used throughout the galliform literature. It is also what
# step 05 reads if you run it with CHR_NAMING=homology.
#
# WHY IT IS NEEDED
#   bLagMut1 is a Sanger/VGP curated assembly. Its chromosomes are named
#   SUPER_1 ... SUPER_38, SUPER_Z, SUPER_W, and that numbering is a LENGTH
#   RANKING, not homology to chicken — the autosome numbers run in near-perfect
#   descending size order and step over the Z, which is the signature. Step 03
#   checks for this and says so.
#
#   So ptarmigan SUPER_6 is not chicken chr6, and renaming it chr_6 would make
#   "chr_6" mean one chromosome in the chicken-scaffolded tree and a different
#   one in the ptarmigan-scaffolded tree. Nothing would fail; every cross-tree
#   comparison would simply be wrong. This job computes the correspondence so
#   the names can be made to agree.
#
# WHAT IT ANSWERS DIRECTLY
#   The open question from the chicken run is whether grouse carry chicken chr6
#   and chr8 as a single chromosome. The suggestive arithmetic is:
#       chicken chr_6 + chr_8   65.80 Mb   (measured in GALGAL_GRCg7b_ref)
#       grouse single chromosome 65.85 Mb  (median across 46 haplotypes)
#       ptarmigan SUPER_5        65.88 Mb
#   Three numbers within 0.08 Mb. That is a hypothesis, not a result — two
#   chromosomes can share a length without sharing ancestry. If this job reports
#       FUSION: ptarmigan 5 = chicken 6 + 8
#   then the fusion is shared with another Tetraoninae genome assembled by
#   someone else, which is independent of anything our pipeline did.
#
#   Watch also for chicken chr2 (~149 Mb): ptarmigan's largest autosome after
#   SUPER_1 is 110 Mb, so chr2 very likely comes back as a SPLIT. Whether that
#   is real or an artefact of either assembly is a separate question.
#
# WHAT IT DOES NOT ANSWER
#   A correspondence says the two genomes differ. It does not say which lineage
#   changed. "Chicken fused 6 and 8" and "ptarmigan split them" produce the same
#   table; distinguishing them needs an outgroup. The literature on GGA4 (a known
#   fusion, retained as two chromosomes in turkey) is the kind of evidence that
#   does the polarising.
#
# USAGE
#   sbatch 04_chromosome_homology.sh
#
#   Then read qc/chromosome_homology/ptarmigan_vs_chicken.summary.tsv, and if
#   the correspondence looks sound, you may run step 05 with CHR_NAMING=homology (its
#   default) so the grouse assemblies inherit chicken-comparable names.
#
# RUNTIME: one minimap2 pass over two ~1 Gb genomes. Under an hour at 24 cpus.
# =============================================================================
#SBATCH --job-name=chr_homology
#SBATCH --output=logs/%x_%j.out
#SBATCH --error=logs/%x_%j.err
#SBATCH -A fnrdewoody
#SBATCH -t 08:00:00
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=24
# 64G: minimap2 asm20 on a ~1 Gb target. The index dominates; the Python step
# that follows holds only per-chromosome sums and is negligible.
#SBATCH --mem=64G
#SBATCH -p cpu
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=blackan@purdue.edu

set -euo pipefail

module unload anaconda 2>/dev/null || true
ml biocontainers
ml minimap2
ml samtools/1.22.1

# unset LD_PRELOAD: RCAC's XALT library is injected this way and fails on some
# nodes (GLIBC mismatch), which can kill subshells under `set -e`.
unset LD_PRELOAD

# =============================================================================
# USER-DEFINED VARIABLES
# =============================================================================
PROJECT_DIR="${CLUSTER_SCRATCH}/GROUSE/grouse_asm"
REF_DIR="${PROJECT_DIR}/ref"
OUT_DIR="${PROJECT_DIR}/qc/chromosome_homology"

# The chicken reference, downloaded by step 01. This is the TARGET: the
# karyotype whose chromosome numbers the galliform literature uses.
CHICKEN_ASM_DIR="GCF_016699485.2_bGalGal1.mat.broiler.GRCg7b"
CHICKEN_FASTA="${REF_DIR}/${CHICKEN_ASM_DIR}_genomic.fna"
CHICKEN_CHR_MAP="${REF_DIR}/${CHICKEN_ASM_DIR}.chr_map.tsv"

# The ptarmigan reference, resolved by step 03. This is the QUERY.
ASM_DIR_RECORD="${REF_DIR}/ptarmigan_asm_dir.txt"

SCRIPT_DIR="${SLURM_SUBMIT_DIR:-$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)}"
HOMOLOGY_SCRIPT="${SCRIPT_DIR}/04_chromosome_homology.py"

# ---- Thresholds -------------------------------------------------------------
# A chicken chromosome counts as "substantial" for a ptarmigan chromosome when
# it holds at least this fraction of that chromosome's matched bases. 0.15 is
# deliberately permissive: a genuine fusion of a 36 Mb and a 30 Mb chromosome
# splits roughly 55/45, while a shared repeat family contributes a percent or
# two, so anything in between is worth seeing rather than hiding.
MIN_FRAC=0.15

# Ignore short alignment blocks. Repeat-driven hits are what scatter alignments
# across non-homologous chromosomes, and they are overwhelmingly short.
MIN_BLOCK=5000
MIN_MAPQ=10

# asm20 as everywhere else in this project: these are between-genus comparisons.
MINIMAP_PRESET="asm20"

THREADS=$SLURM_CPUS_PER_TASK

mkdir -p logs "$OUT_DIR"

# =============================================================================
# PRE-FLIGHT
# =============================================================================
if [[ ! -f "$HOMOLOGY_SCRIPT" ]]; then
    echo "ERROR: companion script not found: ${HOMOLOGY_SCRIPT}"
    echo "       04_chromosome_homology.py must sit next to this file."
    exit 1
fi
for BIN in minimap2 samtools python3; do
    if ! command -v "$BIN" >/dev/null 2>&1; then
        echo "ERROR: '${BIN}' not on PATH."
        exit 1
    fi
done

if [[ ! -s "$ASM_DIR_RECORD" ]]; then
    echo "ERROR: ${ASM_DIR_RECORD} not found."
    echo "       Run 03_ptarmigan_reference.sh first."
    exit 1
fi
ASM_DIR="$(head -n1 "$ASM_DIR_RECORD")"
PTARMIGAN_FASTA="${REF_DIR}/${ASM_DIR}_genomic.fna"
PTARMIGAN_CHR_MAP="${REF_DIR}/${ASM_DIR}.chr_map.tsv"

for F in "$PTARMIGAN_FASTA" "$PTARMIGAN_CHR_MAP"; do
    if [[ ! -s "$F" ]]; then
        echo "ERROR: missing ${F} — run 03_ptarmigan_reference.sh first."
        exit 1
    fi
done
for F in "$CHICKEN_FASTA" "$CHICKEN_CHR_MAP"; do
    if [[ ! -s "$F" ]]; then
        echo "ERROR: missing ${F}"
        echo "       The chicken reference and its chromosome map come from"
        echo "       01_download_chicken_reference.sh."
        exit 1
    fi
done

for F in "$PTARMIGAN_FASTA" "$CHICKEN_FASTA"; do
    if [[ ! -s "${F}.fai" || "${F}.fai" -ot "$F" ]]; then
        echo ">>> Indexing $(basename "$F")"
        samtools faidx "$F"
    fi
done

echo ">>> query  (ptarmigan): ${ASM_DIR}"
echo ">>> target (chicken)  : ${CHICKEN_ASM_DIR}"
echo ">>> ptarmigan chromosomes: $(wc -l < "$PTARMIGAN_CHR_MAP")"
echo ">>> chicken   chromosomes: $(wc -l < "$CHICKEN_CHR_MAP")"

# =============================================================================
# WHOLE-GENOME ALIGNMENT
# =============================================================================
PAF="${OUT_DIR}/ptarmigan_vs_chicken.paf"

if [[ -s "$PAF" && "$PAF" -nt "$PTARMIGAN_FASTA" && "$PAF" -nt "$CHICKEN_FASTA" ]]; then
    echo ">>> alignment already exists and is current — skipping"
else
    echo ">>> minimap2 -x ${MINIMAP_PRESET} (target=chicken, query=ptarmigan)"
    echo "    started: $(date)"
    # Argument order is target then query, which is what makes PAF columns 1-4
    # the PTARMIGAN sequence and columns 6-9 the CHICKEN sequence. The Python
    # step depends on that orientation; swapping it silently inverts every call.
    minimap2 -x "$MINIMAP_PRESET" -t "$THREADS" --secondary=no \
        "$CHICKEN_FASTA" "$PTARMIGAN_FASTA" > "${PAF}.part"
    mv "${PAF}.part" "$PAF"
    echo "    finished: $(date)"
fi
echo ">>> PAF: ${PAF} ($(wc -l < "$PAF") records)"

# =============================================================================
# CORRESPONDENCE
# =============================================================================
PY="$(command -v python3)"
if ! "$PY" -c "import sys; sys.exit(0)" >/dev/null 2>&1; then
    echo "ERROR: python3 is not usable."
    exit 1
fi

echo ""
"$PY" "$HOMOLOGY_SCRIPT" \
    --paf            "$PAF" \
    --query-chr-map  "$PTARMIGAN_CHR_MAP" \
    --target-chr-map "$CHICKEN_CHR_MAP" \
    --query-fai      "${PTARMIGAN_FASTA}.fai" \
    --target-fai     "${CHICKEN_FASTA}.fai" \
    --out-prefix     "${OUT_DIR}/ptarmigan_vs_chicken" \
    --query-label    "ptarmigan" \
    --target-label   "chicken" \
    --min-frac       "$MIN_FRAC" \
    --min-block      "$MIN_BLOCK" \
    --min-mapq       "$MIN_MAPQ"

echo ""
echo "============================================================"
echo ">>> Complete: $(date)"
echo "  summary        : ${OUT_DIR}/ptarmigan_vs_chicken.summary.tsv"
echo "  full pairs     : ${OUT_DIR}/ptarmigan_vs_chicken.correspondence.tsv"
echo "  rename map     : ${OUT_DIR}/ptarmigan_vs_chicken.rename_map.tsv"
echo ""
echo "  READ THE SUMMARY BEFORE RUNNING STEP 09."
echo "  Step 05 defaults to CHR_NAMING=native. Set CHR_NAMING=homology there to"
echo "  above, so that chr_N means the same chromosome in both trees. Set"
echo "  consume the rename map below instead of the reference's own"
echo "  size-rank numbering instead (not comparable to the chicken tree)."
echo "============================================================"
