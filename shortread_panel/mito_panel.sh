#!/bin/bash
# =============================================================================
# mito_panel.sh -- mitochondrial lineage (LEPC vs GRPC/STGR) of every bird in the hybrid-index run,
# by running mito_read_check.sh over the panel in parallel (Gautschi).
#   LEPC carries its own maternal lineage; GRPC and STGR share another, so frac_LEPC ~1 = LEPC mother
#   line, ~0 = GRPC (or STGR) mother line. With the hybrid index this gives the maternal species of F1s
#   and the direction of backcrossing.
#
#   bash mito_panel.sh prep      # login node: diagnostic sites (first bird), ID list
#   bash mito_panel.sh submit    # array over NTASK slices -> combine
#   bash mito_panel.sh status
# Output: ${OUT}/mito_panel.tsv (one row per bird)
# Birds: ${HI_OUT}/ids.txt from hybrid_index.sh setup (the 452 birds, 2026 hold-backs excluded).
# =============================================================================
#SBATCH -A dewoody
#SBATCH -p cpu
set -euo pipefail
unset LD_PRELOAD
REPO_DIR="${SLURM_SUBMIT_DIR:-$PWD}"                       # run from shortread_panel/
SCRIPT="${REPO_DIR}/mito_panel.sh"
HI_OUT="${HI_OUT:-${CLUSTER_SCRATCH}/GROUSE/hybrid_index}"
export OUT="${OUT:-${CLUSTER_SCRATCH}/GROUSE/grouse_asm/mitogenomes/read_check_panel}"
NTASK="${NTASK:-24}"
QOS="${QOS:-normal}"; WALL="${WALL:-08:00:00}"
die() { echo "ERROR: $*" >&2; exit 1; }

case "${1:-}" in
  prep)
    mkdir -p "${OUT}/tables" "${OUT}/logs"
    [[ -s "${HI_OUT}/ids.txt" ]] || die "no ${HI_OUT}/ids.txt (run hybrid_index.sh setup)"
    cp "${HI_OUT}/ids.txt" "${OUT}/ids.txt"
    first=$(head -1 "${OUT}/ids.txt")
    TABLE="${OUT}/tables/prep.tsv" bash "${REPO_DIR}/mito_read_check.sh" "$first"
    echo "diagnostic sites: $(($(wc -l < "${OUT}/diagnostic_sites.tsv") - 1)); birds: $(wc -l < "${OUT}/ids.txt")" ;;
  task)
    i=${SLURM_ARRAY_TASK_ID:?}
    mapfile -t ids < <(awk -v n="$NTASK" -v i="$i" '(NR - 1) % n == i' "${OUT}/ids.txt")
    TABLE="${OUT}/tables/part_${i}.tsv" bash "${REPO_DIR}/mito_read_check.sh" "${ids[@]}" ;;
  combine)
    awk 'FNR == 1 && NR != 1 {next} {print}' "${OUT}"/tables/part_*.tsv > "${OUT}/mito_panel.tsv"
    echo "birds: $(($(wc -l < "${OUT}/mito_panel.tsv") - 1)) -> ${OUT}/mito_panel.tsv" ;;
  submit)
    [[ -s "${OUT}/diagnostic_sites.tsv" ]] || die "run: bash $0 prep"
    J=$(sbatch --parsable -A dewoody -p cpu -q "$QOS" -t "$WALL" --array=0-$((NTASK - 1)) -c 8 --mem=16G \
        --export=ALL -J mito_panel -o "${OUT}/logs/%x_%A_%a.out" "$SCRIPT" task)
    J2=$(sbatch --parsable -A dewoody -p cpu -q "$QOS" -t 00:10:00 -c 1 --mem=2G --dependency=afterany:"$J" \
        --export=ALL -J mito_combine -o "${OUT}/logs/%x_%j.out" "$SCRIPT" combine)
    echo "array $J ($NTASK tasks) -> combine $J2" ;;
  status)
    echo "birds done: $(cat "${OUT}"/tables/part_*.tsv 2>/dev/null | grep -vc '^ID' || true) / $(wc -l < "${OUT}/ids.txt")"
    squeue -u "$USER" -n mito_panel,mito_combine 2>/dev/null | head ;;
  *) sed -n '2,14p' "$0" | sed 's/^# \{0,1\}//'; exit 1 ;;
esac
