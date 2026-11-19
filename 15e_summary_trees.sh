#!/bin/bash
# =============================================================================
# STEP 13e: ONE SUMMARY TREE FOR THE AUTOSOMES AND ONE FOR THE Z (report Figure 17)
#
# Run after 13_thex_phylogeny.sh (and 13d, if used) on Negishi. Nothing is realigned:
# the step-13 window trees (work/iqtree/<chrom>/<window>.contree) are summarised.
#
#   collect   concatenate the window trees: autosomes (every chromosome except Z) and Z
#   astral    wASTRAL (ASTER; Zhang & Mirarab 2022): coalescent summary tree that weights
#             each window tree by its UFBoot support and branch lengths, so windows whose
#             relationships are unresolved count for little. Two runs per partition:
#               haplotype level  46 haplotypes + ptarmigan as tips (local posterior support)
#               species level    tips mapped to LEPC / GRPC / STGR / Lmuta (-a), with the
#                                quartet frequencies of the three possible species trees
#   gcf       gene concordance factors (IQ-TREE --gcf): for every branch of the haplotype
#             tree, the percentage of window trees that contain it
#   submit    all of the above as one standby job (~1 h)
#   setup     once: conda env with ASTER at ${PROJECT_DIR}/envs/aster_env
#
# Why not concatenation or a majority-rule consensus: a concatenated ML tree of 1 Gb gives
# 100% bootstrap to every node whether or not windows agree, and can be wrong when lineage
# sorting is incomplete; a majority-rule consensus of 9,000 discordant window trees collapses
# to a polytomy among the prairie-chickens. ASTRAL is consistent under incomplete lineage
# sorting, and gCF shows how much of the genome supports each branch.
#
# Outputs: work/summary_trees/
#   {autosomes,chrZ}.windows.trees        input window trees (one Newick per line)
#   {autosomes,chrZ}.astral.tre           haplotype-level tree (local PP on nodes; internal
#                                         branch lengths in coalescent units)
#   {autosomes,chrZ}.species.astral.tre   species-level tree, full annotation (q1 = share of
#                                         quartets supporting the shown branch)
#   {autosomes,chrZ}.gcf.cf.tree          haplotype tree with gCF on nodes  -> 13f plot
#   {autosomes,chrZ}.gcf.cf.stat          gCF / gDF1 / gDF2 / gDFP per branch
#
# USAGE (Negishi, repo root)
#   bash 13e_summary_trees.sh setup
#   bash 13e_summary_trees.sh submit
#   then copy work/summary_trees/*.cf.tree and the 12d mito tree to the laptop and run
#   Rscript 13f_plot_summary_trees.R autosomes.gcf.cf.tree chrZ.gcf.cf.tree mito_tree.contree
# =============================================================================
#SBATCH -A fnrdewoody
#SBATCH -p cpu
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=blackan@purdue.edu

set -euo pipefail
unset LD_PRELOAD

REPO_DIR="${SLURM_SUBMIT_DIR:-$PWD}"
SCRIPT="${REPO_DIR}/13e_summary_trees.sh"
PROJECT_DIR="${PROJECT_DIR:-${CLUSTER_SCRATCH}/GROUSE/grouse_thex}"
WORK="${PROJECT_DIR}/work"
IQ_DIR="${WORK}/iqtree"
OUT="${WORK}/summary_trees"
THEX_ENV="${THEX_ENV:-${PROJECT_DIR}/envs/thex_env}"     # IQ-TREE 2 for --gcf (from step 13 setup)
ASTER_ENV="${ASTER_ENV:-${PROJECT_DIR}/envs/aster_env}"
OUTGROUP_NAME="${OUTGROUP_NAME:-Lmuta}"
Z_DIR_RE="${Z_DIR_RE:-^chr[-_]Z$}"                       # step-13 chromosome folder of the Z
THREADS="${THREADS:-${SLURM_CPUS_PER_TASK:-8}}"
ACCOUNT="${ACCOUNT:-fnrdewoody}"
PARTITION="${PARTITION:-cpu}"
CPUS="${CPUS:-32}"

MODE="${1:-}"
die() { echo "ERROR: $*" >&2; exit 1; }
log() { echo "[$(date '+%F %T')] $*"; }

find_astral() {
    export PATH="${ASTER_ENV}/bin:${PATH}"
    ASTRAL_BIN="$(command -v wastral || command -v astral4 || command -v astral || true)"
    [[ -n "$ASTRAL_BIN" ]] || die "no wastral/astral in ${ASTER_ENV}; run: bash $0 setup"
    ASTRAL_HELP="$("$ASTRAL_BIN" -h 2>&1 || true)"
    # wASTRAL weights each window tree by its UFBoot values (0-100) and branch lengths
    # (hybrid weighting, --mode 1). -S = bootstrap support scale; NOT -B, which is abayes [0,1].
    SUPPORT_FLAG=()
    if [[ "$(basename "$ASTRAL_BIN")" == wastral && "$ASTRAL_HELP" == *"--bootstrap"* ]]; then
        SUPPORT_FLAG=(-S)
    else
        log "NOTE: $(basename "$ASTRAL_BIN") is not wASTRAL; window trees are not weighted by support"
    fi
}

stage_setup() {
    module load conda >/dev/null 2>&1 || module load anaconda >/dev/null 2>&1 || true
    command -v conda >/dev/null || die "no conda: try 'module avail conda anaconda'"
    export CONDA_PKGS_DIRS="${PROJECT_DIR}/envs/pkgs"
    if [[ ! -d "$ASTER_ENV" ]]; then
        local solver; solver="$(command -v mamba || command -v conda)"
        "$solver" create -y -p "$ASTER_ENV" --override-channels -c conda-forge -c bioconda aster
    fi
    find_astral
    log "astral: $ASTRAL_BIN  support flag: ${SUPPORT_FLAG[*]:-none}"
    [[ -x "${THEX_ENV}/bin/iqtree2" ]] || die "IQ-TREE 2 not in ${THEX_ENV} (run step 13 setup)"
}

stage_collect() {
    [[ -d "$IQ_DIR" ]] || die "no step-13 trees in $IQ_DIR"
    mkdir -p "$OUT"
    : > "${OUT}/autosomes.windows.trees"; : > "${OUT}/chrZ.windows.trees"
    : > "${OUT}/autosomes.windows.list"; : > "${OUT}/chrZ.windows.list"
    local d c part
    for d in "$IQ_DIR"/*/; do
        c="$(basename "$d")"
        if [[ "$c" =~ $Z_DIR_RE ]]; then part=chrZ; else part=autosomes; fi
        find "$d" -maxdepth 1 -name '*.contree' | sort -V | while read -r f; do
            tr -d '\n' < "$f" >> "${OUT}/${part}.windows.trees"; echo >> "${OUT}/${part}.windows.trees"
            echo "$f" >> "${OUT}/${part}.windows.list"
        done
    done
    log "window trees: autosomes $(wc -l < "${OUT}/autosomes.windows.trees"), Z $(wc -l < "${OUT}/chrZ.windows.trees")"
    log "(report Table 3: 8,851 autosomal and 686 Z windows)"
    # tip -> species map for the species-level runs (LEPC_F1234_hap1 -> LEPC)
    cat "${OUT}"/*.windows.trees | tr '(),' '\n\n\n' | sed 's/:.*//' | grep -E '^[A-Za-z]' | sort -u \
        | awk -v og="$OUTGROUP_NAME" '{ s = ($1 == og) ? og : $1; sub(/_.*/, "", s); print $1"\t"s }' \
        > "${OUT}/species_map.tsv"
    log "tips: $(wc -l < "${OUT}/species_map.tsv") ($(cut -f2 "${OUT}/species_map.tsv" | sort | uniq -c | tr -s ' ' | tr '\n' ' '))"
}

stage_astral() {
    find_astral
    local p
    for p in autosomes chrZ; do
        [[ -s "${OUT}/${p}.windows.trees" ]] || die "run collect first"
        log "$p: haplotype-level tree"
        "$ASTRAL_BIN" ${SUPPORT_FLAG[@]+"${SUPPORT_FLAG[@]}"} -t "$THREADS" --root "$OUTGROUP_NAME" -u 1 -i "${OUT}/${p}.windows.trees" \
            -o "${OUT}/${p}.astral.tre" 2> "${OUT}/${p}.astral.log"
        log "$p: species-level tree"
        "$ASTRAL_BIN" ${SUPPORT_FLAG[@]+"${SUPPORT_FLAG[@]}"} -t "$THREADS" --root "$OUTGROUP_NAME" -u 2 -a "${OUT}/species_map.tsv" \
            -i "${OUT}/${p}.windows.trees" -o "${OUT}/${p}.species.astral.tre" 2> "${OUT}/${p}.species.astral.log"
        echo "  $(cat "${OUT}/${p}.species.astral.tre")"
    done
}

stage_gcf() {
    local IQ="${THEX_ENV}/bin/iqtree2" p
    for p in autosomes chrZ; do
        # topology only: IQ-TREE does not need ASTRAL's support labels or coalescent lengths
        sed -E 's/:[^,();]+//g; s/\)[^,();]+/)/g' "${OUT}/${p}.astral.tre" > "${OUT}/${p}.astral.topology.tre"
        "$IQ" -t "${OUT}/${p}.astral.topology.tre" --gcf "${OUT}/${p}.windows.trees" \
            --prefix "${OUT}/${p}.gcf" -redo --quiet
        log "$p: gCF -> ${OUT}/${p}.gcf.cf.tree"
    done
}

stage_submit() {
    mkdir -p "${REPO_DIR}/logs"
    sbatch -A "$ACCOUNT" -p "$PARTITION" -q standby -t 04:00:00 --export=ALL -c "$CPUS" \
        -J thex_summary_trees -o "${REPO_DIR}/logs/%x_%j.out" "$SCRIPT" all
}

case "$MODE" in
    setup)   stage_setup ;;
    collect) stage_collect ;;
    astral)  stage_astral ;;
    gcf)     stage_gcf ;;
    all)     stage_collect; stage_astral; stage_gcf
             log "done; copy ${OUT}/*.cf.tree and *.species.astral.tre to the laptop" ;;
    submit)  stage_submit ;;
    *) sed -n '2,42p' "$0" | sed 's/^# \{0,1\}//'; exit 1 ;;
esac
