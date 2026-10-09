#!/bin/bash
# =============================================================================
# Relabel chr_W in MALE (ZZ) haplotypes as unplaced sequence — both trees.
#
# Males have no W chromosome. Anything RagTag placed on the reference W in a
# male is repeat / Z-gametolog sequence, so it is renamed scaffold_N (next free
# number, never colliding with an existing name). Females are never touched.
#
#   >= 50 kb  -> stays in <prefix>.pseudo_chr.fasta as scaffold_N
#   <  50 kb  -> moved to <prefix>.unplaced_short.fasta (same rule as step 09)
#
# Also renames chr_W in the matching .liftoff.gff3, saves the original chr_W
# sequence to <prefix>.chrW_original.fa, and logs every change to
# <tree>/W_relabel_log.tsv.
#
# By default it processes EVERY final-assembly directory it finds among
# CANDIDATE_TREES below (chicken_guided/final, plus the ptarmigan tree wherever
# it currently lives) and reports any it cannot find.
#
# Safe to rerun: a haplotype with no chr_W left is skipped.
#
# USAGE (login node is fine; it takes minutes):
#   MALES="F5457" bash relabel_male_W.sh                       # test, both trees
#   bash relabel_male_W.sh                                     # all males, both trees
#   TREE=chicken_guided/final bash relabel_male_W.sh           # one tree only
# =============================================================================
set -euo pipefail

PROJECT_DIR="${CLUSTER_SCRATCH}/GROUSE/grouse_asm"
MIN="${MIN:-50000}"

# Final-assembly directories, relative to PROJECT_DIR. Every one that exists is
# processed. Several ptarmigan locations are listed so this keeps working if
# that tree is moved too. TREE=<path> overrides the whole list.
CANDIDATE_TREES=(
    "chicken_guided/final"
    "final_ptarmigan"
    "ptarmigan_guided/final_ptarmigan"
    "ptarmigan_guided/final"
)
if [[ -n "${TREE:-}" ]]; then CANDIDATE_TREES=("$TREE"); fi

MALES="${MALES:-F5540 F5541 F5542 F5543 F5544 F5545 F5546 F5457 F5462 F5463 F5468 F5472 F5478 F5480 F5485 F5502 F5597}"
# The six ZW females. Hard stop if any is in MALES by mistake.
FEMALES="F5595 F5596 F5598 F5599 F5600 F5503"

cd "$PROJECT_DIR"

# Load samtools inside the script, the same way the pipeline's sbatch scripts
# do. A module loaded in the login shell does not always reach a child
# `bash script.sh` (RCAC biocontainer commands can be shell aliases/functions,
# which a child shell does not inherit), so do not rely on the caller's
# environment. Lmod's own scripts reference unset variables, so relax -u
# around the load.
if ! command -v samtools >/dev/null 2>&1; then
    set +u
    module unload anaconda >/dev/null 2>&1 || true
    module load biocontainers  >/dev/null 2>&1 || true
    module load samtools/1.22.1 >/dev/null 2>&1 || true
    set -u
fi
# RCAC's XALT tracker is injected via LD_PRELOAD and can kill subshells on some
# nodes; the pipeline scripts drop it for the same reason.
unset LD_PRELOAD

if ! command -v samtools >/dev/null 2>&1; then
    echo "ERROR: samtools still not available after 'module load biocontainers samtools/1.22.1'."
    echo "       What the shell sees:"
    type -a samtools 2>&1 | sed 's/^/         /' || true
    echo "       Is the module function available here?  type module -> $(type -t module || echo none)"
    echo "       Workaround: run it as a job instead —  sbatch --wrap='bash relabel_male_W.sh' -A dewoody -p cpu -t 1:00:00"
    exit 1
fi
echo ">>> samtools: $(command -v samtools)"
for S in $MALES; do
    for X in $FEMALES; do
        [[ "$S" == "$X" ]] && { echo "ERROR: ${S} is female — remove it from MALES"; exit 1; }
    done
done

# -----------------------------------------------------------------------------
relabel_tree() {
    local T="$1"
    local LOG="${T}/W_relabel_log.tsv"
    [[ -f "$LOG" ]] || printf 'haplotype\told\tnew\tlength_bp\tdestination\tdate\n' > "$LOG"

    local N_DONE=0 S H F L B SHORT N NEW DEST N_BEFORE N_AFTER EXPECT GFF
    for S in $MALES; do
      for H in hap1 hap2; do
        F=$(ls "${T}"/*_"${S}"_"${H}".pseudo_chr.fasta 2>/dev/null | head -1 || true)
        if [[ -z "$F" || ! -s "$F" ]]; then echo "  skip ${S} ${H}: no FASTA in ${T}"; continue; fi
        [[ -s "${F}.fai" && "${F}.fai" -nt "$F" ]] || samtools faidx "$F"

        L=$(awk '$1=="chr_W"{print $2}' "${F}.fai")
        if [[ -z "$L" ]]; then echo "  skip ${S} ${H}: no chr_W"; continue; fi

        B="${F%.pseudo_chr.fasta}"
        SHORT="${B}.unplaced_short.fasta"
        [[ -f "$SHORT" ]] || : > "$SHORT"

        # next free scaffold number across BOTH the main and short files
        N=$(grep -ho '^>scaffold_[0-9]*' "$F" "$SHORT" 2>/dev/null | sed 's/^>scaffold_//' | sort -n | tail -1 || true)
        NEW="scaffold_$(( ${N:-0} + 1 ))"
        if grep -q "^>${NEW}\$" "$F" "$SHORT" 2>/dev/null; then
            echo "ERROR: ${NEW} already exists in ${B} — stopping"; exit 1
        fi

        # keep the original W sequence (small) instead of a full-genome backup
        samtools faidx "$F" chr_W > "${B}.chrW_original.fa"
        sed "1s/.*/>${NEW}/" "${B}.chrW_original.fa" > "${B}.W.tmp.fa"

        awk '$1!="chr_W"{print $1}' "${F}.fai" > "${B}.keep.tmp"
        N_BEFORE=$(wc -l < "${F}.fai")
        samtools faidx "$F" -r "${B}.keep.tmp" > "${F}.tmp"

        if (( L >= MIN )); then
            cat "${B}.W.tmp.fa" >> "${F}.tmp"; DEST="main"; EXPECT=$N_BEFORE
        else
            cat "${B}.W.tmp.fa" >> "$SHORT"; DEST="unplaced_short"; EXPECT=$(( N_BEFORE - 1 ))
        fi

        mv "${F}.tmp" "$F"
        rm -f "${B}.W.tmp.fa" "${B}.keep.tmp" "${F}.fai"
        samtools faidx "$F"
        if [[ -s "$SHORT" ]]; then rm -f "${SHORT}.fai"; samtools faidx "$SHORT"; fi

        # sanity: no chr_W left, sequence count reconciles
        if grep -q '^>chr_W$' "$F"; then echo "ERROR: chr_W still present in ${F}"; exit 1; fi
        N_AFTER=$(wc -l < "${F}.fai")
        if (( N_AFTER != EXPECT )); then
            echo "ERROR: ${F}: ${N_BEFORE} sequences before, ${N_AFTER} after (expected ${EXPECT})"; exit 1
        fi

        GFF="${B}.liftoff.gff3"
        if [[ -s "$GFF" ]]; then
            awk -F'\t' -v OFS='\t' -v n="$NEW" '$1=="chr_W"{$1=n}1' "$GFF" > "${GFF}.tmp" && mv "${GFF}.tmp" "$GFF"
        fi

        printf '%s\tchr_W\t%s\t%s\t%s\t%s\n' "$(basename "$B")" "$NEW" "$L" "$DEST" "$(date +%F)" >> "$LOG"
        printf '  %-22s chr_W (%9s bp) -> %-14s [%s]\n' "$(basename "$B")" "$L" "$NEW" "$DEST"
        N_DONE=$(( N_DONE + 1 ))
      done
    done

    echo ""
    echo "  relabelled ${N_DONE} haplotype(s). Log: ${PROJECT_DIR}/${LOG}"
    echo "  haplotypes still carrying chr_W (should be females only):"
    grep -l '^>chr_W$' "${T}"/*_hap[12].pseudo_chr.fasta 2>/dev/null \
        | xargs -r -n1 basename | sed 's/\.pseudo_chr\.fasta$//; s/^/      /'
}
# -----------------------------------------------------------------------------

N_TREES=0
for T in "${CANDIDATE_TREES[@]}"; do
    if ls "${T}"/*_hap[12].pseudo_chr.fasta >/dev/null 2>&1; then
        echo ""
        echo "============================================================"
        echo ">>> ${PROJECT_DIR}/${T}"
        echo "============================================================"
        relabel_tree "$T"
        N_TREES=$(( N_TREES + 1 ))
    else
        echo ">>> not found (skipped): ${T}"
    fi
done

echo ""
if (( N_TREES == 0 )); then
    echo "ERROR: no final-assembly directory found. Pass one with TREE=<path relative to ${PROJECT_DIR}>"
    exit 1
fi
echo ">>> done: ${N_TREES} tree(s) processed"
