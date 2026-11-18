#!/bin/bash
# =============================================================================
# NCBI SUBMISSION: .sqn files with per-gap linkage evidence
# Run after ncbi_submission_prep.sh (uses its files/, rename_maps/ and
# haplotype_table.tsv). NCBI accepts N-runs as gaps automatically only when all
# gaps share one linkage evidence, and no longer accepts AGP files, so each
# assembly is converted to .sqn with table2asn from a feature table that gives
# every gap its own evidence (see ncbi_gap_tables.py):
#   YaHS gaps (Hi-C scaffolding)                -> proximity-ligation
#   RagTag joins (ordered against ptarmigan)    -> align-xgenus
#   RagTag joins supported by the Hi-C test     -> align-xgenus + proximity-ligation
#
# STAGES (run from accessory/)
#   sbatch --export=ALL,STAGE=tables ncbi_make_sqn.sh       # .fsa + .tbl per assembly, gap_summary.tsv
#   ONLY=F5540 sbatch --export=ALL,STAGE=tables ncbi_make_sqn.sh   # one-bird test
#   bash ncbi_make_sqn.sh sqn                                # array: table2asn on all 46 (needs TEMPLATE)
#   bash ncbi_make_sqn.sh check                              # errors in the .val/.stats reports
#
# TEMPLATE: a submission template (.sbt) made at
#   https://submit.ncbi.nlm.nih.gov/genbank/template/submission/
#   (contact and author details), saved on Gautschi, e.g. ${SUB}/template.sbt
# =============================================================================
#SBATCH --job-name=ncbi_sqn
#SBATCH -A dewoody
#SBATCH -p cpu
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --mem=32G
#SBATCH -t 08:00:00
#SBATCH --output=logs/%x_%A_%a.out
#SBATCH --error=logs/%x_%A_%a.err
set -euo pipefail
unset LD_PRELOAD

PROJECT_DIR="${CLUSTER_SCRATCH}/GROUSE/grouse_asm"
SUB="${SUB:-${PROJECT_DIR}/ncbi_submission}"
SQN="${SQN:-${SUB}/sqn}"
RAGTAG_DIR="${RAGTAG_DIR:-${PROJECT_DIR}/ragtag_ptarmigan}"
JOINS_DIR="${JOINS_DIR:-${PROJECT_DIR}/qc_ptarmigan/join_support}"
TEMPLATE="${TEMPLATE:-${SUB}/template.sbt}"
BIN="${PROJECT_DIR}/bin"
ONLY="${ONLY:-}"
STAGE="${STAGE:-${1:-}}"
HERE="${SLURM_SUBMIT_DIR:-$(pwd)}"
mkdir -p logs "$SQN" "$BIN"

get_table2asn() {
    T2A="$(command -v table2asn || true)"
    if [[ -z "$T2A" ]]; then
        T2A="${BIN}/table2asn"
        if [[ ! -x "$T2A" ]]; then
            curl -sSL https://ftp.ncbi.nlm.nih.gov/asn1-converters/by_program/table2asn/linux64.table2asn.gz \
                | gunzip > "$T2A" && chmod +x "$T2A"
        fi
    fi
    "$T2A" -version | head -1
}

case "$STAGE" in
tables)
    python3 "${HERE}/ncbi_gap_tables.py" --sub "$SUB" --ragtag-dir "$RAGTAG_DIR" --joins-dir "$JOINS_DIR" \
        --out "$SQN" --only "$ONLY"
    column -t -s $'\t' "${SQN}/gap_summary.tsv" | cut -c1-200
    ;;
sqn)
    if [[ -z "${SLURM_ARRAY_TASK_ID:-}" ]]; then       # login node: check inputs, submit the array
        [[ -s "$TEMPLATE" ]] || { echo "ERROR: no template at ${TEMPLATE} (see header)"; exit 1; }
        get_table2asn
        ls -d "${SQN}"/F*_hap? > "${SQN}/dirs.txt"
        N=$(wc -l < "${SQN}/dirs.txt")
        sbatch --export=ALL,STAGE=sqn --array=0-$((N - 1))%12 "$0"
        echo "Submitted ${N} table2asn tasks"
        exit 0
    fi
    get_table2asn
    D=$(sed -n "$((SLURM_ARRAY_TASK_ID + 1))p" "${SQN}/dirs.txt")
    echo ">>> ${D}  $(date)"
    # -M n: genome checks; -Z: discrepancy report; gaps come from the .tbl assembly_gap features
    "$T2A" -indir "$D" -outdir "$D" -t "$TEMPLATE" -M n -Z -V v
    ls -lh "$D"/*.sqn
    echo ">>> done $(date)"
    ;;
check)
    set +e +o pipefail                                  # report, never stop, if files are missing
    n=0
    for d in "${SQN}"/F*_hap?; do
        [[ -d "$d" ]] || continue
        n=$((n + 1)); b=$(basename "$d")
        nsqn=$(ls "$d"/*.sqn 2>/dev/null | wc -l)
        nval=$(ls "$d"/*.val 2>/dev/null | wc -l)
        nerr=$(cat "$d"/*.val 2>/dev/null | grep -c -E '^(ERROR|REJECT)')
        printf "%s\tsqn=%s\tval=%s\terrors=%s\n" "$b" "$nsqn" "$nval" "$nerr"
    done
    echo "assemblies checked: ${n} (in ${SQN})"
    echo; echo "Error types across all files:"
    cat "${SQN}"/F*_hap?/*.val 2>/dev/null | grep -E '^(ERROR|REJECT)' | sed -E 's/^([A-Z]+): valid \[([^]]+)\].*/\1 \2/' | sort | uniq -c | sort -rn | head -20
    ;;
*)
    echo "Usage: sbatch --export=ALL,STAGE=tables $0 | bash $0 sqn | bash $0 check" >&2; exit 1 ;;
esac
