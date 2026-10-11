#!/bin/bash
# =============================================================================
# rerun_without_ONT.sh -- re-assemble F5545 (LEPC) and F5597 (GRPC) from HiFi +
# Hi-C only (no ONT ultra-long reads), like the other 21 birds, and carry the
# new assemblies through every downstream step that the report depends on.
#
# Why: with ONT-UL (hifiasm --ul) these four haplotypes had Merqury QV 51.7-55.9
# against 62.2-64.7 for the HiFi-only haplotypes, with 6-15x more assembly-only
# k-mers, while contiguity (N50, scaffold count) was no better.
#
# Run on a Gautschi LOGIN node, in three stages:
#
#   bash rerun_without_ONT.sh plan      # read-only: indices, files to archive, jobs
#   bash rerun_without_ONT.sh prepare   # archive old outputs (mv, nothing deleted)
#   bash rerun_without_ONT.sh submit    # submit the dependency chain below
#
# Chain (assembly repo = prairie_grouse, pangenome repo = grouse_pangenome):
#   A  02_assembly_scaffold_array.sh   array (2 tasks)  hifiasm -> yahs -> .hic   ~10 d
#   B  05_ragtag_liftoff_array.sh      array (2 tasks)  RagTag/Liftoff/QUAST/Hi-C ~3 d
#   C  11_join_support.sh              array (2 tasks)  Hi-C support of each join
#   D  05b_relabel_male_W.sh           2 males          chr_W -> scaffold_N
#   E  07_busco_array.sh               array (2 tasks)  BUSCO + Z/W/MT split
#   F  10_sex_chromosome_check.sh      all 46 (fast)    Table S3
#   G  09_busco_plots.sh               all (~4 h)       Figure S6 (barplot)
#   H  01_break_unsupported_joins.sh   all 46 (< 1 h)   pangenome inputs
#   I  02_assembly_qc_array.sh         array (2 tasks)  Merqury + HiFi depth
# then by hand: 02b_collate_qc.sh, 02c_qc_report_table.sh, local plots.
#
# BEFORE RUNNING
#   1. git pull in BOTH repos on the cluster: both assembly_manifest.tsv files
#      must already have ont_ul = NA for F5545 and F5597 (this script checks).
#   2. If you made manual Juicebox (.review.assembly) edits to F5545 or F5597,
#      they are archived, not reapplied: redo them on the new .hic/.assembly
#      before step B (submit with  STAGES=B  afterwards; see below).
#   3. Set ASM_REPO / PG_REPO below to the cluster clones of the two repos.
# =============================================================================
set -euo pipefail

SAMPLES=(F5545 F5597)
PROJECT_DIR="${CLUSTER_SCRATCH}/GROUSE/grouse_asm"
ASM_REPO="${ASM_REPO:-${HOME}/prairie_grouse}"       # cluster clone, assembly repo
PG_REPO="${PG_REPO:-${HOME}/grouse_pangenome}"       # cluster clone, pangenome repo
ACCOUNT="${ACCOUNT:-fnrdewoody}"                     # for the small wrapper jobs
ARCHIVE="${PROJECT_DIR}/archive_with_ONT"            # old outputs are moved here
STAGES="${STAGES:-A}"                                # first stage to submit (A..I)

# Top-level directories that never hold per-sample ASSEMBLY outputs:
# raw reads, adapter-filtered HiFi (reusable: does not depend on ONT), references,
# environments, tools, earlier archives, and the chicken-era trees.
KEEP_DIRS=(raw hifi_filtered ref conda_envs tools archive_with_ONT chicken_guided)

die() { echo "ERROR: $*" >&2; exit 1; }
mode="${1:-plan}"

# ---- array index of each sample, exactly as the array scripts compute it ----
index_of() {   # $1 manifest, $2 sample
    local i=0 line
    while IFS= read -r line; do
        [[ "${line%%$'\t'*}" == "$2" ]] && { echo "$i"; return; }
        i=$((i+1))
    done < <(grep -v '^#' "$1" | tail -n +2 | grep -v '^[[:space:]]*$')
    die "$2 not in $1"
}
check_manifest() {  # ont_ul (column 6) must be NA for both samples
    local m="$1" s v
    [[ -f "$m" ]] || die "manifest not found: $m"
    for s in "${SAMPLES[@]}"; do
        v=$(awk -F'\t' -v s="$s" '$1==s {print $6}' "$m" | tr -d '\r')
        [[ "$v" == "NA" ]] || die "$m: ont_ul for $s is '$v' (must be NA; git pull?)"
    done
}
ASM_MAN="${ASM_REPO}/assembly_manifest.tsv"; PG_MAN="${PG_REPO}/assembly_manifest.tsv"
check_manifest "$ASM_MAN"; check_manifest "$PG_MAN"
IDX=(); for s in "${SAMPLES[@]}"; do IDX+=("$(index_of "$ASM_MAN" "$s")"); done
PIDX=(); for s in "${SAMPLES[@]}"; do PIDX+=("$(index_of "$PG_MAN" "$s")"); done
ARR=$(IFS=,; echo "${IDX[*]}"); PARR=$(IFS=,; echo "${PIDX[*]}")

# ---- everything to archive: any path under PROJECT_DIR named for the sample --
# A matching directory is moved whole (-prune). In pangenome_qc/<SAMPLE>/ the
# HiFi read FASTA and meryl database are KEPT (they depend only on the reads,
# and meryl counting is the slow part); everything else there is archived.
targets() {
    local s="$1" d keep=()
    for d in "${KEEP_DIRS[@]}"; do keep+=(-path "${PROJECT_DIR}/${d}" -o); done
    find "$PROJECT_DIR" -mindepth 1 \( "${keep[@]}" -path "${PROJECT_DIR}/pangenome_qc/${s}" \) -prune \
         -o -name "*${s}*" -print -prune
    if [[ -d "${PROJECT_DIR}/pangenome_qc/${s}" ]]; then
        find "${PROJECT_DIR}/pangenome_qc/${s}" -mindepth 1 -maxdepth 1 \
             ! -name "${s}.hifi.fa*" ! -name "*.meryl" -print
    fi
}

case "$mode" in
plan)
    echo "Samples: ${SAMPLES[*]}"
    echo "Assembly-repo array indices: ${ARR}    pangenome-repo array indices: ${PARR}"
    echo "Manifests OK (ont_ul = NA): $ASM_MAN, $PG_MAN"
    echo; echo "Would move to ${ARCHIVE}/ (same relative layout):"
    for s in "${SAMPLES[@]}"; do targets "$s" | sed "s|^${PROJECT_DIR}/|  |"; done
    echo; echo "Next: bash $0 prepare   then   bash $0 submit"
    ;;
prepare)
    mkdir -p "$ARCHIVE"
    n=0
    for s in "${SAMPLES[@]}"; do
        while IFS= read -r p; do
            rel="${p#${PROJECT_DIR}/}"; dest="${ARCHIVE}/${rel}"
            mkdir -p "$(dirname "$dest")"
            [[ -e "$dest" ]] && die "already archived: $dest (move it away first)"
            mv "$p" "$dest"; n=$((n+1))
        done < <(targets "$s")
    done
    echo "Archived ${n} paths to ${ARCHIVE}"
    echo "Paths still matching (expected: none; raw/, hifi_filtered/ and pangenome_qc reads/meryl are kept by design):"
    for s in "${SAMPLES[@]}"; do targets "$s"; done | sed 's/^/  /'
    ;;
submit)
    sb() { sbatch --parsable "$@"; }
    dep() { [[ -n "${1:-}" ]] && echo "--dependency=afterok:$1" || true; }
    order=(A B C D E F G H I); start=0
    for i in "${!order[@]}"; do [[ "${order[$i]}" == "$STAGES" ]] && start=$i; done
    run() { local st="$1"; for i in "${!order[@]}"; do [[ "${order[$i]}" == "$st" ]] && (( i >= start )) && return 0; done; return 1; }
    A=; B=; C=; D=; E=; F=; G=; H=; I=
    cd "$ASM_REPO"; mkdir -p logs
    run A && A=$(sb --array="$ARR" 02_assembly_scaffold_array.sh)
    run B && B=$(sb $(dep "$A") --array="$ARR" --export=ALL,CHR_NAMING=native,DO_HIC=true 05_ragtag_liftoff_array.sh)
    run C && C=$(sb $(dep "$B") --array="$ARR" 11_join_support.sh)
    run D && D=$(sb $(dep "$C") -A "$ACCOUNT" -p cpu -t 01:00:00 -c 2 --mem=8G -J relabel_W \
                   -o logs/relabel_W_%j.out --wrap "cd '$ASM_REPO' && MALES='${SAMPLES[*]}' TREE=final_ptarmigan bash 05b_relabel_male_W.sh")
    run E && E=$(sb $(dep "$D") --array="$ARR" 07_busco_array.sh)
    run F && F=$(sb $(dep "$E") 10_sex_chromosome_check.sh)
    run G && G=$(sb $(dep "$E") 09_busco_plots.sh)
    cd "$PG_REPO"; mkdir -p logs
    run H && H=$(sb $(dep "$D") 01_break_unsupported_joins.sh)
    run I && I=$(sb $(dep "$H") --array="$PARR" 02_assembly_qc_array.sh)
    printf "Submitted (stage: job id)\n  A 02 assembly: %s\n  B 05 ragtag:   %s\n  C 11 joins:    %s\n  D 05b W:       %s\n  E 07 BUSCO:    %s\n  F 10 sex:      %s\n  G 09 plots:    %s\n  H pg01 inputs: %s\n  I pg02 QC:     %s\n" \
        "${A:--}" "${B:--}" "${C:--}" "${D:--}" "${E:--}" "${F:--}" "${G:--}" "${H:--}" "${I:--}"
    cat <<MSG

When I has finished (from ${PG_REPO}):
  bash 02b_collate_qc.sh
  bash 02c_qc_report_table.sh          # new pangenome_qc_summary.tsv -> 02d_plot_qc_report.R
New QUAST values for Figure S5 (scaffolds, total length, N50 of the final assemblies):
  for s in ${SAMPLES[*]}; do for h in hap1 hap2; do
    f=\$(ls -d ${PROJECT_DIR}/qc_ptarmigan/quast/*\${s}*\${h}*/report.tsv 2>/dev/null | head -1)
    echo "\$s \$h  \$f"; grep -E '^(# contigs\t|Total length\t|N50\t|L50\t|GC)' "\$f"; done; done
MSG
    ;;
*) die "usage: bash $0 plan|prepare|submit" ;;
esac
