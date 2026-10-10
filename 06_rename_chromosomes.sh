#!/bin/bash
# =============================================================================
# RENAME THE CHROMOSOMES IN THE PTARMIGAN-SCAFFOLDED ASSEMBLIES
# Step 06 — run after 05, and after 07/08/11 if you have already run them.
#
# ONLY NEEDED if step 05 was run with CHR_NAMING=homology. Its default is
# CHR_NAMING=native, which produces these names directly and makes this step a
# no-op. It exists because the delivered assemblies were built the other way
# round, and because it is the safe way to rename a tree that already has
# derived files (GFF3, AGP, CRAM, BUSCO tables) pointing at the old names.
#
# WHY
#   The chr_ names currently carry a comparison rather than naming a chromosome:
#     chr_6_8    "spans chicken 6 and 8"
#     chr_2a/2b  "pieces of chicken 2"
#     chr_4a/4b  "pieces of chicken 4"      (4b is GGA4p)
#     chr_u30    "could not be matched to chicken"
#   Every one of those is a statement about CHICKEN. They were the right names
#   while choosing between references. They are the wrong names on a finished
#   assembly, where an identifier should be short, unique and stable, and the
#   homology belongs in a table beside it.
#
# SCHEMES (set SCHEME below)
#   native   chr_1 ... chr_38, chr_Z, chr_W, chr_MT — the ptarmigan reference's
#            own chromosome numbering. Every chromosome gets a plain name, no
#            chr_u*, no compound names, no a/b. Numbering is fixed by the single
#            reference, so it is IDENTICAL across all 46 haplotypes — which is
#            the property that matters and that per-assembly size-ranking would
#            not give you. This is the default and the recommended scheme.
#            Caveat: chr_N here does NOT mean chicken chr_N. The homology table
#            this script writes is what carries that, and it must travel with
#            the assemblies.
#
#   swap-ab  Keep the chicken-homology names but letter the split segments by
#            POSITION along the chicken chromosome instead of by size, so chr_4a
#            is GGA4p rather than the larger piece. Minimal change; use this if
#            you want to stay with chicken-homologous names.
#
#   custom   Use the TSV at CUSTOM_MAP (old<TAB>new) as-is.
#
# SAFETY
#   Dry run by default — set APPLY=true to write. The mapping is checked for
#   collisions before anything is touched, renaming is simultaneous (so swaps
#   work), every file is written to a .part and moved on success, and a reverse
#   map is written so the whole thing can be undone.
#
# NOT HANDLED
#   Pretext maps (.pretext) are binary and cannot be edited; regenerate them if
#   you want the new names — the command is printed at the end. QUAST reports
#   keep the old names; they are descriptive output, not inputs to anything.
#
# USAGE
#   bash 06_rename_chromosomes.sh              # dry run, prints the plan
#   APPLY=true bash 06_rename_chromosomes.sh   # do it
# =============================================================================
#SBATCH --job-name=grouse_rename_chr
#SBATCH --output=logs/%x_%j.out
#SBATCH --error=logs/%x_%j.err
#SBATCH -A fnrdewoody
#SBATCH -t 04:00:00
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --mem=16G
#SBATCH -p cpu
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=blackan@purdue.edu

set -euo pipefail

module unload anaconda 2>/dev/null || true
ml biocontainers
ml samtools/1.22.1
unset LD_PRELOAD

# =============================================================================
# USER-DEFINED VARIABLES
# =============================================================================
PROJECT_DIR="${CLUSTER_SCRATCH}/GROUSE/grouse_asm"
REF_DIR="${PROJECT_DIR}/ref"
REF_TAG="ptarmigan"

SCHEME="${SCHEME:-native}"          # native | swap-ab | custom
APPLY="${APPLY:-false}"             # true actually writes
CUSTOM_MAP="${CUSTOM_MAP:-}"        # only for SCHEME=custom

# Rewrite the Hi-C CRAM headers too. Sequence data is untouched, so the CRAM's
# per-sequence M5 checksums stay valid — only the @SQ SN names move. Set false
# to skip (the CRAMs are only needed for contact maps and step 11).
DO_CRAM="${DO_CRAM:-true}"

HOMOLOGY_DIR="${PROJECT_DIR}/qc/chromosome_homology"
HOMOLOGY_MAP="${HOMOLOGY_DIR}/ptarmigan_vs_chicken.rename_map.tsv"
ASM_DIR_RECORD="${REF_DIR}/ptarmigan_asm_dir.txt"

OUT_DIR="${PROJECT_DIR}/qc_${REF_TAG}/chromosome_naming"
SCRIPT_DIR="${SLURM_SUBMIT_DIR:-$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)}"
RENAME_PY="${SCRIPT_DIR}/06_rename_chromosomes.py"

mkdir -p logs "$OUT_DIR"

# =============================================================================
# PRE-FLIGHT
# =============================================================================
[[ -f "$RENAME_PY" ]] || { echo "ERROR: 06_rename_chromosomes.py must sit next to this file."; exit 1; }
for BIN in samtools python3 awk; do
    command -v "$BIN" >/dev/null 2>&1 || { echo "ERROR: '${BIN}' not on PATH."; exit 1; }
done
[[ -s "$ASM_DIR_RECORD" ]] || { echo "ERROR: run 08 first (${ASM_DIR_RECORD} missing)."; exit 1; }
ASM_DIR="$(head -n1 "$ASM_DIR_RECORD")"
PTARM_CHR_MAP="${REF_DIR}/${ASM_DIR}.chr_map.tsv"
[[ -s "$PTARM_CHR_MAP" ]] || { echo "ERROR: ${PTARM_CHR_MAP} missing."; exit 1; }
[[ -s "$HOMOLOGY_MAP" ]]  || { echo "ERROR: ${HOMOLOGY_MAP} missing — run 08b."; exit 1; }

MAP="${OUT_DIR}/rename.${SCHEME}.tsv"
TABLE="${OUT_DIR}/chromosome_homology_table.tsv"

# =============================================================================
# BUILD THE MAPPING
# Both inputs are keyed on the ptarmigan RefSeq accession, so the current name
# and the native name are joined through it rather than guessed.
#   ${PTARM_CHR_MAP}  accession -> native number   (NC_064437.1 -> 5)
#   ${HOMOLOGY_MAP}   accession -> chicken label   (NC_064437.1 -> 6_8)
# =============================================================================
echo ">>> scheme : ${SCHEME}"
echo ">>> ref    : ${ASM_DIR}"

case "$SCHEME" in
native)
    awk -F'\t' '
        FNR==NR { if ($0 !~ /^#/ && $1 != "") nat[$1] = $2; next }
        /^#/ || $1 == "" { next }
        {
            acc = $1; lab = $2
            if (!(acc in nat)) next
            old = "chr_" lab
            new = "chr_" nat[acc]
            print old "\t" new
        }
    ' "$PTARM_CHR_MAP" "$HOMOLOGY_MAP" | sort -u > "${MAP}.part"
    mv "${MAP}.part" "$MAP"
    ;;
swap-ab)
    # Position-ordered lettering for split segments only; everything else keeps
    # its current name. Derived from the homology summary, which records which
    # segment sits where along the chicken chromosome.
    SUM="${HOMOLOGY_DIR}/ptarmigan_vs_chicken.summary.tsv"
    [[ -s "$SUM" ]] || { echo "ERROR: ${SUM} missing — rerun 08b."; exit 1; }
    echo "  deriving a/b swaps from ${SUM}"
    awk -F'\t' 'NR>1 && $3=="split" {print $6}' "$SUM" | sort > "${MAP}.labels"
    if [[ ! -s "${MAP}.labels" ]]; then
        echo "ERROR: no split rows in ${SUM}; nothing to swap."
        exit 1
    fi
    # Pair up a<->b within each target chromosome.
    awk '{
        lab=$0; sub(/^chr_/,"",lab)
        base=lab; suf=substr(lab,length(lab),1)
        sub(/[a-z]$/,"",base)
        if (suf ~ /[a-z]/) grp[base]=grp[base] " " suf
        full[base "_" suf]=$0
    }
    END {
        for (b in grp) {
            n=split(grp[b],s," ")
            for (i=1;i<=n;i++) {
                j=n-i+1
                print full[b "_" s[i]] "\t" "chr_" b s[j]
            }
        }
    }' "${MAP}.labels" | sort -u > "${MAP}.part"
    rm -f "${MAP}.labels"
    mv "${MAP}.part" "$MAP"
    ;;
custom)
    [[ -s "$CUSTOM_MAP" ]] || { echo "ERROR: SCHEME=custom needs CUSTOM_MAP set to a TSV."; exit 1; }
    cp "$CUSTOM_MAP" "$MAP"
    ;;
*)
    echo "ERROR: SCHEME must be native, swap-ab or custom (got '${SCHEME}')"
    exit 1
    ;;
esac

if [[ ! -s "$MAP" ]]; then
    echo "ERROR: mapping came out empty — check ${PTARM_CHR_MAP} and ${HOMOLOGY_MAP}"
    exit 1
fi
echo ">>> mapping: ${MAP} ($(grep -vc '^#' "$MAP") entries)"

# =============================================================================
# THE HOMOLOGY TABLE — what has to travel with the assemblies
# Once the names no longer encode chicken homology, this file is the only place
# that records it. Write it whether or not the rename is applied.
# =============================================================================
awk -F'\t' -v OFS='\t' '
    FNR==NR { if ($0 !~ /^#/ && $1 != "") nat[$1]=$2; next }
    /^#/ || $1=="" { next }
    FNR==1 && NR!=1 { }
    { acc=$1; lab=$2; if (!(acc in nat)) next
      print "chr_" nat[acc], "chr_" lab, acc, nat[acc], lab }
' "$PTARM_CHR_MAP" "$HOMOLOGY_MAP" \
 | sort -V > "${TABLE}.part"
{
  printf 'assembly_chr\tprevious_name\tptarmigan_accession\tptarmigan_chr\tchicken_homolog\n'
  cat "${TABLE}.part"
} > "$TABLE"
rm -f "${TABLE}.part"
echo ">>> homology table: ${TABLE}"
echo "    KEEP THIS WITH THE ASSEMBLIES — once the names are native, it is the"
echo "    only record of which chicken chromosome each one corresponds to."

# =============================================================================
# APPLY
# =============================================================================
PY_ARGS=(--project-dir "$PROJECT_DIR" --ref-tag "$REF_TAG" --map "$MAP")
if [[ "$APPLY" == true ]]; then
    PY_ARGS+=(--apply --reverse-map-out "${OUT_DIR}/rename.${SCHEME}.reverse.tsv")

    # Say up front how much I/O this is. Every FASTA is read and rewritten in
    # full, because header lines change length and cannot be patched in place.
    FA_BYTES=$(du -cb "${PROJECT_DIR}/final_${REF_TAG}"/*.pseudo_chr.fasta \
                      "${PROJECT_DIR}/final_autosomes_${REF_TAG}"/*.autosome_chr.fasta \
                      2>/dev/null | tail -n1 | cut -f1 || echo 0)
    if [[ "${FA_BYTES:-0}" -gt 0 ]]; then
        echo ""
        echo ">>> about to rewrite $(awk -v b="$FA_BYTES" 'BEGIN{printf "%.1f", b/1e9}') GB of FASTA"
        echo "    (read + write, so roughly double that in I/O). Expect minutes,"
        echo "    not seconds. Progress is printed per file."
        echo "    Free space here: $(df -h "${PROJECT_DIR}" | awk 'NR==2{print $4}')"
    fi

    # A .part left by an interrupted run is harmless — the original was never
    # replaced — but say so rather than leaving it to be discovered later.
    STALE=$(find "${PROJECT_DIR}/final_${REF_TAG}" "${PROJECT_DIR}/final_autosomes_${REF_TAG}" \
                 -maxdepth 1 -name '*.part' 2>/dev/null | wc -l)
    if [[ "${STALE:-0}" -gt 0 ]]; then
        echo "    note: ${STALE} leftover .part file(s) from an earlier interrupted run;"
        echo "          the originals are intact and these will be overwritten."
    fi
fi

echo ""
python3 "$RENAME_PY" "${PY_ARGS[@]}"

# =============================================================================
# FASTA INDICES AND CRAM HEADERS — the parts that need samtools
# =============================================================================
if [[ "$APPLY" == true ]]; then
    echo ""
    echo ">>> re-indexing FASTAs"
    N=0
    for F in "${PROJECT_DIR}/final_${REF_TAG}"/*.pseudo_chr.fasta \
             "${PROJECT_DIR}/final_autosomes_${REF_TAG}"/*.autosome_chr.fasta; do
        [[ -s "$F" ]] || continue
        rm -f "${F}.fai"
        samtools faidx "$F"          # also fails loudly on any duplicate name
        N=$((N+1))
    done
    echo "    ${N} FASTA(s) re-indexed"

    if [[ "$DO_CRAM" == true ]]; then
        echo ""
        echo ">>> rewriting Hi-C CRAM headers"
        M=0
        for CRAM in "${PROJECT_DIR}/qc_${REF_TAG}"/*.hic2final.cram; do
            [[ -s "$CRAM" ]] || continue
            samtools view -H "$CRAM" > "${CRAM}.hdr"
            awk -v OFS='\t' -v map="$MAP" '
                BEGIN { while ((getline l < map) > 0) {
                            if (l ~ /^#/) continue
                            split(l, a, "\t"); if (a[1] != "") m[a[1]] = a[2] } }
                /^@SQ/ { for (i=1;i<=NF;i++) if ($i ~ /^SN:/) {
                             n = substr($i,4); if (n in m) $i = "SN:" m[n] } }
                { print }
            ' "${CRAM}.hdr" > "${CRAM}.hdr.new"
            # Sequence data is untouched, so the M5 checksums remain correct.
            samtools reheader "${CRAM}.hdr.new" "$CRAM" > "${CRAM}.part"
            mv "${CRAM}.part" "$CRAM"
            samtools index "$CRAM"
            rm -f "${CRAM}.hdr" "${CRAM}.hdr.new"
            M=$((M+1))
        done
        echo "    ${M} CRAM(s) reheadered and re-indexed"
    fi

    # ---- verification ----
    echo ""
    echo ">>> verifying: no old names left in any assembly FASTA"
    # An old name is only "left over" if it is NOT also somebody's new name.
    # These mappings are full of reshuffles — chr_2a becomes chr_3 while the old
    # chr_3 becomes chr_2 — so chr_3 is still in the file, legitimately, as a
    # different chromosome. Checking old names blindly reports those as failures
    # and makes a correct rename look broken.
    awk -F'\t' '!/^#/ && $2 != "" {print $2}' "$MAP" | sort -u > "${MAP}.newnames"
    LEFT=0
    while IFS=$'\t' read -r OLD NEW; do
        [[ "$OLD" == \#* || -z "${OLD:-}" || "$OLD" == "$NEW" ]] && continue
        if grep -qxF "$OLD" "${MAP}.newnames"; then continue; fi
        if grep -lq "^>${OLD}$" "${PROJECT_DIR}/final_${REF_TAG}"/*.pseudo_chr.fasta 2>/dev/null; then
            echo "    ! '${OLD}' still present and is not a new name — genuine leftover"
            LEFT=$((LEFT+1))
        fi
    done < "$MAP"
    rm -f "${MAP}.newnames"

    # Independent check that nothing was merged: the number of distinct names in
    # each assembly must be unchanged. samtools faidx above would already have
    # failed on a duplicate, so this is belt and braces.
    for F in "${PROJECT_DIR}/final_${REF_TAG}"/*.pseudo_chr.fasta; do
        [[ -s "$F" ]] || continue
        TOTN=$(grep -c '^>' "$F")
        UNIQ=$(grep '^>' "$F" | sort -u | wc -l)
        if (( TOTN != UNIQ )); then
            echo "    ! $(basename "$F"): ${TOTN} sequences but only ${UNIQ} distinct names"
            LEFT=$((LEFT+1))
        fi
    done
    if (( LEFT == 0 )); then
        echo "    clean"
    else
        echo "    ${LEFT} old name(s) remain — investigate before using these assemblies"
        exit 1
    fi

    echo ""
    echo ">>> chromosomes now in a representative assembly:"
    REP=$(ls "${PROJECT_DIR}/final_${REF_TAG}"/*.pseudo_chr.fasta 2>/dev/null | head -n1)
    [[ -n "$REP" ]] && awk -F'\t' '$1 ~ /^chr_/ {printf "      %-10s %10.2f Mb\n", $1, $2/1e6}' "${REP}.fai" | head -50
fi

# =============================================================================
# DONE
# =============================================================================
echo ""
echo "============================================================"
if [[ "$APPLY" == true ]]; then
    echo ">>> Rename applied: $(date)"
    echo "  mapping       : ${MAP}"
    echo "  undo with     : CUSTOM_MAP=${OUT_DIR}/rename.${SCHEME}.reverse.tsv \\"
    echo "                    SCHEME=custom APPLY=true bash 06_rename_chromosomes.sh"
    echo "  homology table: ${TABLE}"
    echo ""
    echo "  Pretext contact maps still carry the OLD names and cannot be edited."
    echo "  To regenerate one:"
    echo "    samtools view -h -T <final.fasta> <qc_${REF_TAG}/S.H.hic2final.cram> \\"
    echo "      | PretextMap -o <qc_${REF_TAG}/S.H.pretext> --sortby length --sortorder descend --mapq 20"
    echo ""
    echo "  If 10/11 already ran, rerun 09_busco_plots.sh so the"
    echo "  karyotype and synteny figures carry the new names."
else
    echo ">>> DRY RUN — nothing was changed: $(date)"
    echo "  mapping that WOULD be applied: ${MAP}"
    echo "  homology table written anyway: ${TABLE}"
    echo ""
    echo "  Review the mapping, then:"
    echo "    APPLY=true bash 06_rename_chromosomes.sh"
fi
echo "============================================================"
