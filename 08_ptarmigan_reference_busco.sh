#!/bin/bash
# =============================================================================
# SLURM JOB: ROCK PTARMIGAN REFERENCE AS A CONTRAST ASSEMBLY
# Step 08 — requires
# 03_ptarmigan_reference.sh, and is intended to run before
# 09_busco_plots.sh.
#
# NOT an array job. One genome, two steps.
#
# The point of this script is to make the ptarmigan reference look like just
# another assembly, so steps 10 and 11 pick it up with no special-casing:
#
#   1. Rename  — write ${FINAL_DIR}/LAGMUT_bLagMut1_ref.pseudo_chr.fasta, the
#                ptarmigan reference with its RefSeq accessions replaced by the
#                SAME chr_* names the grouse assemblies now carry, via the same
#                NCBI assembly-report map step 03 built and step 05 renamed
#                against. Unplaced sequences become scaffold_N under the same
#                minimum-length filter, so the contrast is like-for-like.
#   2. BUSCO   — run into ${QC_DIR}/busco/LAGMUT_bLagMut1_ref/, the same tree as
#                the ptarmigan-scaffolded grouse runs, same lineage.
#
# Once this finishes, run 09_busco_plots.sh and the ptarmigan appears
# automatically: its own completeness bar, its own karyotype plot, and a synteny
# panel against one representative of each grouse species.
#
# WHY THIS ONE MATTERS MORE THAN THE CHICKEN EQUIVALENT DID
#   The chicken contrast answered "how do our assemblies compare against the
#   standard galliform reference". This one answers the question that raised:
#   the ptarmigan karyotype plot drawn beside the grouse karyotype plots is a
#   direct test of whether chr_6/chr_8 and the chr_4 ~19 Mb join are
#   Tetraoninae-wide features or artefacts of our assemblies —
#     - ptarmigan chr_6 ~65 Mb with no separate chr_8 -> the fusion is ancestral
#       to Tetraoninae, our assemblies are right, chicken is the outlier
#     - ptarmigan chr_6 ~36 Mb plus chr_8 ~30 Mb      -> ptarmigan shares
#       chicken's state, so the fusion is specific to Tympanuchus: still a real
#       finding, but a different one
#   Note this is the reference's own karyotype, i.e. independent evidence — it
#   does not depend on anything our pipeline did.
#
# The renamed FASTA also works as input to 07_busco_array.sh's
# splitting logic if you ever want the reference's chr_W/chr_Z/chr_MT pulled out.
#
# USAGE:
#   sbatch 08_ptarmigan_reference_busco.sh
#
# WHY RENAMING MATTERS: without it the ptarmigan sequences are NC_*/NW_
# accessions while the grouse assemblies are chr_1 etc. The synteny plot joins
# BUSCOs on id but labels chromosomes by name, and the karyotype builder keeps
# only chr_*-prefixed sequences — so an unrenamed reference would produce an
# empty karyotype and an unreadable synteny panel.
# =============================================================================
#SBATCH --job-name=ptarmigan_ref_busco
#SBATCH --output=logs/%x_%j.out
#SBATCH --error=logs/%x_%j.err
#SBATCH -A fnrdewoody
#SBATCH -t 1-00:00:00
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=24
# 64G: same budget as the grouse BUSCO runs — the ptarmigan genome is the same
# size class (~1.0 Gb), and metaeuk dominates the footprint.
#SBATCH --mem=64G
#SBATCH -p cpu
#SBATCH --mail-type=BEGIN,END,FAIL
#SBATCH --mail-user=blackan@purdue.edu

# =============================================================================
# ENVIRONMENT SETUP
# =============================================================================
set -euo pipefail

module unload anaconda 2>/dev/null || true

ml biocontainers
ml busco
ml samtools/1.22.1

# unset LD_PRELOAD: RCAC's XALT usage-tracking library is injected via
# LD_PRELOAD and fails on some nodes (GLIBC_2.33/2.34 mismatch), which can
# kill subshells under `set -e`. It's accounting only — safe to drop.
unset LD_PRELOAD

# =============================================================================
# USER-DEFINED VARIABLES
# =============================================================================
PROJECT_DIR="${CLUSTER_SCRATCH}/GROUSE/grouse_asm"
REF_DIR="${PROJECT_DIR}/ref"

# Must match REF_TAG in 05_ragtag_liftoff_array.sh.
REF_TAG="ptarmigan"

FINAL_DIR="${PROJECT_DIR}/final_${REF_TAG}"
QC_DIR="${PROJECT_DIR}/qc_${REF_TAG}"
BUSCO_DIR="${QC_DIR}/busco"

# Resolved by step 03 rather than hardcoded — the assembly directory name is the
# one thing about this accession that cannot be known without asking NCBI, and
# having two scripts guess it independently is how they come to disagree.
ASM_DIR_RECORD="${REF_DIR}/ptarmigan_asm_dir.txt"

# The label the ptarmigan carries through every downstream plot. Must match
# <SPECIES>_<SAMPLE>_<HAP> with HAP = 'ref', which is what 09_busco_plots.py
# accepts for an unphased reference genome.
#
# IMPORTANT: the middle field must contain NO underscore. 09_busco_plots.py
# parses directory names with ([A-Z]+)_([A-Za-z0-9.]+)_(hap\d+|ref), so a label
# built straight from the NCBI assembly name (e.g. "bLagMut1_primary") would
# fail to parse and the reference would be silently dropped from every plot.
# That is why this is a fixed short label rather than derived from ASM_DIR.
REF_LABEL="LAGMUT_bLagMut1_ref"

# Same threshold as step 05, so the reference's unplaced pile is filtered on the
# same terms as the grouse assemblies. chr_* sequences are always kept
# regardless of length (galliform microchromosomes are legitimately small).
MIN_UNPLACED_SCAFFOLD_LEN=50000

BUSCO_LINEAGE="aves_odb10"
BUSCO_MODE="genome"
BUSCO_DOWNLOAD_PATH="${PROJECT_DIR}/busco_downloads"

THREADS=$SLURM_CPUS_PER_TASK

mkdir -p logs "$FINAL_DIR" "$BUSCO_DIR" "$BUSCO_DOWNLOAD_PATH"

# Same guarded helper as 07_busco_array.sh: on a first run the BUSCO
# run directory does not exist, `find` exits 1, pipefail carries it through the
# pipe, and under `set -e` a failing command substitution in an assignment kills
# the script before BUSCO is ever called.
find_busco_summary() {
    local dir="$1"
    [[ -d "$dir" ]] || { printf ''; return 0; }
    find "$dir" -maxdepth 1 -name 'short_summary.specific.*.txt' -print 2>/dev/null \
        | head -n1 || true
}

for BIN in busco samtools awk; do
    if ! command -v "$BIN" >/dev/null 2>&1; then
        echo "ERROR: '${BIN}' not on PATH after module load."
        exit 1
    fi
done
echo ">>> busco    : $(command -v busco)  [$(busco --version 2>&1 | head -n1 || true)]"
echo ">>> samtools : $(command -v samtools)"

# =============================================================================
# RESOLVE THE REFERENCE
# Everything here is produced by step 03. This script deliberately does NOT
# re-download or rebuild any of it: if the map used to rename the grouse
# assemblies in step 05 and the map used to rename the reference here were ever
# built from different inputs, the chr_* names would stop meaning the same thing
# and every synteny panel would be quietly wrong.
# =============================================================================
if [[ ! -s "$ASM_DIR_RECORD" ]]; then
    echo "ERROR: ${ASM_DIR_RECORD} not found."
    echo "       Run 03_ptarmigan_reference.sh first."
    exit 1
fi
REF_ASM_DIR="$(head -n1 "$ASM_DIR_RECORD")"
if [[ -z "$REF_ASM_DIR" ]]; then
    echo "ERROR: ${ASM_DIR_RECORD} is empty."
    exit 1
fi

SRC_FASTA="${REF_DIR}/${REF_ASM_DIR}_genomic.fna"
REF_REPORT="${REF_DIR}/${REF_ASM_DIR}_assembly_report.txt"
REF_CHR_MAP="${REF_DIR}/${REF_ASM_DIR}.chr_map.tsv"

echo ">>> Reference: ${REF_ASM_DIR}"

if [[ ! -s "$SRC_FASTA" ]]; then
    echo "ERROR: ptarmigan reference FASTA not found: ${SRC_FASTA}"
    echo "       Run 03_ptarmigan_reference.sh first."
    exit 1
fi
if [[ ! -s "$REF_CHR_MAP" ]]; then
    echo "ERROR: chromosome map not found: ${REF_CHR_MAP}"
    echo "       Run 03_ptarmigan_reference.sh first."
    exit 1
fi
echo ">>> Chromosome map: $(wc -l < "$REF_CHR_MAP") assembled molecules"

# Catch a label that step 09's parser would silently reject, before spending a
# BUSCO run on a directory name that will be skipped.
if [[ ! "$REF_LABEL" =~ ^[A-Z]+_[A-Za-z0-9.]+_ref$ ]]; then
    echo "ERROR: REF_LABEL='${REF_LABEL}' does not match <SPECIES>_<SAMPLE>_ref"
    echo "       with SAMPLE free of underscores. 09_busco_plots.py would not"
    echo "       parse it and the reference would be dropped from every plot."
    exit 1
fi

# =============================================================================
# STEP 1 — Rename the ptarmigan reference to chr_* / scaffold_N
# =============================================================================
REF_FASTA="${FINAL_DIR}/${REF_LABEL}.pseudo_chr.fasta"
SHORT_FASTA="${FINAL_DIR}/${REF_LABEL}.unplaced_short.fasta"
RENAME_MAP="${FINAL_DIR}/${REF_LABEL}.rename_map.tsv"

echo ""
echo ">>> [1/2] Renaming ptarmigan sequences to chr_* / scaffold_N"

if [[ -s "$REF_FASTA" && -s "${REF_FASTA}.fai" && -f "$SHORT_FASTA" \
      && "$REF_FASTA" -nt "$SRC_FASTA" ]] \
   && ! grep -q '^>N[CW]_' "$REF_FASTA"; then
    echo "  already renamed — skipping"
else
    if [[ ! -s "${SRC_FASTA}.fai" || "${SRC_FASTA}.fai" -ot "$SRC_FASTA" ]]; then
        samtools faidx "$SRC_FASTA"
    fi

    # Accession -> chr_<name> where the assembly report lists it as an
    # assembled molecule; everything else becomes scaffold_N, numbered by
    # descending length so the numbering is stable and meaningful.
    awk -F'\t' '{print $1"\t"$2}' "${SRC_FASTA}.fai" \
        | sort -k2,2 -nr \
        | awk -F'\t' -v chrmap="$REF_CHR_MAP" '
            BEGIN {
                while ((getline line < chrmap) > 0) {
                    split(line, a, "\t")
                    chrname[a[1]] = a[2]
                }
            }
            {
                if ($1 in chrname) print $1"\tchr_"chrname[$1]
                else { scafn++; print $1"\tscaffold_"scafn }
            }
        ' > "$RENAME_MAP"
    echo "  rename map: ${RENAME_MAP} ($(grep -c 'chr_' "$RENAME_MAP") chromosomes)"

    awk -v map="$RENAME_MAP" '
        BEGIN {
            while ((getline line < map) > 0) {
                split(line, a, "\t")
                newname[a[1]] = a[2]
            }
        }
        /^>/ {
            split(substr($0, 2), parts, " ")
            name = parts[1]
            print ">" (name in newname ? newname[name] : name)
            next
        }
        { print }
    ' "$SRC_FASTA" > "${REF_FASTA}.prefilter"

    # Same minimum-length rule as the grouse assemblies: chr_* always kept,
    # short unplaced scaffolds moved aside rather than discarded.
    samtools faidx "${REF_FASTA}.prefilter"
    awk -F'\t' -v min="$MIN_UNPLACED_SCAFFOLD_LEN" \
        '$1 ~ /^chr_/ || $2 >= min {print $1}' "${REF_FASTA}.prefilter.fai" > "${REF_FASTA}.keep.txt"
    awk -F'\t' -v min="$MIN_UNPLACED_SCAFFOLD_LEN" \
        '$1 !~ /^chr_/ && $2 < min {print $1}' "${REF_FASTA}.prefilter.fai" > "${REF_FASTA}.short.txt"

    if [[ ! -s "${REF_FASTA}.keep.txt" ]]; then
        echo "ERROR: length filter would keep nothing — check ${REF_FASTA}.prefilter"
        exit 1
    fi
    samtools faidx "${REF_FASTA}.prefilter" -r "${REF_FASTA}.keep.txt" > "$REF_FASTA"
    if [[ -s "${REF_FASTA}.short.txt" ]]; then
        samtools faidx "${REF_FASTA}.prefilter" -r "${REF_FASTA}.short.txt" > "$SHORT_FASTA"
    else
        : > "$SHORT_FASTA"
    fi
    samtools faidx "$REF_FASTA"

    N_IN=$(wc -l < "${REF_FASTA}.prefilter.fai")
    N_KEEP=$(wc -l < "${REF_FASTA}.keep.txt")
    N_SHORT=$(wc -l < "${REF_FASTA}.short.txt")
    if (( N_KEEP + N_SHORT != N_IN )); then
        echo "ERROR: sequence counts do not reconcile: in ${N_IN}, kept ${N_KEEP}, short ${N_SHORT}"
        exit 1
    fi

    rm -f "${REF_FASTA}.prefilter" "${REF_FASTA}.prefilter.fai"

    echo "  input   : ${N_IN} sequences"
    echo "  kept    : ${N_KEEP} -> $(basename "$REF_FASTA")"
    echo "  short   : ${N_SHORT} -> $(basename "$SHORT_FASTA")"
    echo "  chromosomes present:"
    awk -F'\t' '$1 ~ /^chr_/ {printf "    %-10s %12d bp\n", $1, $2}' "${REF_FASTA}.fai" | head -50
fi

# =============================================================================
# STEP 2 — BUSCO on the renamed ptarmigan reference
# =============================================================================
echo ""
echo ">>> [2/2] BUSCO (${BUSCO_LINEAGE}, ${BUSCO_MODE} mode)"

BUSCO_LINEAGE_DIR="${BUSCO_DOWNLOAD_PATH}/lineages/${BUSCO_LINEAGE}"
if [[ ! -s "${BUSCO_LINEAGE_DIR}/dataset.cfg" ]]; then
    echo ">>> Downloading BUSCO lineage ${BUSCO_LINEAGE}"
    busco --download "$BUSCO_LINEAGE" --download_path "$BUSCO_DOWNLOAD_PATH" \
        || echo "  WARNING: busco --download returned non-zero; checking anyway"
    if [[ ! -s "${BUSCO_LINEAGE_DIR}/dataset.cfg" ]]; then
        echo "ERROR: BUSCO lineage not present at ${BUSCO_LINEAGE_DIR}"
        exit 1
    fi
fi

BUSCO_RUN_DIR="${BUSCO_DIR}/${REF_LABEL}"
BUSCO_SUMMARY=$(find_busco_summary "$BUSCO_RUN_DIR")

if [[ -n "$BUSCO_SUMMARY" && -s "$BUSCO_SUMMARY" && "$BUSCO_SUMMARY" -nt "$REF_FASTA" ]]; then
    echo "  already complete — skipping"
    echo "  $(grep -m1 -E '^[[:space:]]+C:' "$BUSCO_SUMMARY" || true)"
else
    if [[ -d "$BUSCO_RUN_DIR" ]]; then
        echo "  removing incomplete previous run: ${BUSCO_RUN_DIR}"
        rm -rf "$BUSCO_RUN_DIR"
    fi

    busco \
        --in "$REF_FASTA" \
        --out "$REF_LABEL" \
        --out_path "$BUSCO_DIR" \
        --mode "$BUSCO_MODE" \
        --lineage_dataset "$BUSCO_LINEAGE_DIR" \
        --offline \
        --download_path "$BUSCO_DOWNLOAD_PATH" \
        --cpu "$THREADS"

    BUSCO_SUMMARY=$(find_busco_summary "$BUSCO_RUN_DIR")
    if [[ -z "$BUSCO_SUMMARY" || ! -s "$BUSCO_SUMMARY" ]]; then
        echo "ERROR: BUSCO finished without writing a short summary"
        echo "       Inspect ${BUSCO_RUN_DIR}"
        exit 1
    fi
    echo "  $(grep -m1 -E '^[[:space:]]+C:' "$BUSCO_SUMMARY" || true)"
fi

# =============================================================================
# DONE
# =============================================================================
echo ""
echo "============================================================"
echo ">>> Ptarmigan reference ready: $(date)"
echo "  FASTA  : ${REF_FASTA}"
echo "  BUSCO  : ${BUSCO_RUN_DIR}/"
echo ""
echo "  Now rerun the plots to include it:"
echo "    sbatch 09_busco_plots.sh"
echo "============================================================"
