#!/bin/bash
# =============================================================================
# LOCAL ANCESTRY ALONG THE GENOME WITH ANCESTRY_HMM (Corbett-Detig & Nielsen 2017) (Gautschi)
# Unnumbered (shortread_panel/). Builds on hybrid_index.sh: its genotype likelihoods (red/*.npz) and ANGSD
# allele files (gl/*.mafs.gz) give the ancestry-informative sites and reference (LEPC, GRPC) allele frequencies.
#
#   Groups, each run separately (Ancestry_HMM models one admixed population):
#     control  LEPC_ctrl  - 30 Texas/New Mexico LEPC held out of the reference: no GRPC ancestry expected, so
#                           their estimates measure the method's floor (false GRPC ancestry)
#     edge     LEPC_edge  - LEPC 18-71 km from the GRPC range (old, low-level admixture; many short tracts)
#     contact  LEPC_sym   - LEPC from the Kansas contact zone (recent hybrids; long tracts)
#   Each group x chromosome class (A = autosomes, Z) is a separate run; on the Z females are haploid.
#   QC-excluded birds (qc_excluded_birds.tsv) are left out everywhere.
#
#   bash ancestry_hmm.sh setup      # login node: conda env with ancestry_hmm; sites.tsv/.bed
#   bash ancestry_hmm.sh submit     # read counts (array over birds) -> build + run (array over 4 runs) -> summary
#   bash ancestry_hmm.sh status
#   bash ancestry_hmm.sh submit_zfix   # after a full run: Z runs again with the pulse time FIXED at the group's
#                                       # autosomal estimate (runs/<group>_Zfix). A free Z time collapses to a few
#                                       # generations (whole-Z blocks in recent hybrids), which cannot detect short tracts.
#
# Model (per run): two ancestries, LEPC (0) the resident population and GRPC (1) arriving in a single pulse
# whose time is estimated (-p 1 -TIME0 ...); GRPC proportion fixed at the group's mean autosomal hybrid index
# from hybrid_index.tsv. Recombination: RATE_TABLE (per scaffold) or a uniform RATE = 3.11 cM/Mb, the chicken
# sex-averaged genome average (Groenen et al. 2009, Genome Res. 19:510); Z x 2/3 (recombines in males only); Ne NE.
# Output ($OUT): runs/<group>_<class>/<bird>.posterior, ancestry_hmm_windows.tsv (1-Mb windows: mean GRPC
# ancestry across birds), ancestry_hmm_per_bird.tsv, ancestry_hmm.png, and each run's log with the pulse time.
# =============================================================================
#SBATCH -A dewoody
#SBATCH -p cpu
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=blackan@purdue.edu
set -euo pipefail
unset LD_PRELOAD

REPO_DIR="${SLURM_SUBMIT_DIR:-$PWD}"                    # run from shortread_panel/
SCRIPT="${REPO_DIR}/ancestry_hmm.sh"
PY="${REPO_DIR}/ancestry_hmm_prep.py"
HI="${HI:-${CLUSTER_SCRATCH}/GROUSE/hybrid_index}"      # hybrid_index.sh output
OUT="${OUT:-${CLUSTER_SCRATCH}/GROUSE/ancestry_hmm}"
PANEL_GROUPS="${PANEL_GROUPS:-${REPO_DIR}/panel_groups_lepc_grpc.tsv}"
GROUPS_TSV="${OUT}/groups_ahmm.tsv"      # PANEL_GROUPS with CONTROL_N Texas/New Mexico LEPC held out as LEPC_ctrl (setup)
CONTROL_N="${CONTROL_N:-30}"             # held-out unadmixed LEPC: their estimated GRPC ancestry is the method's floor
SEX_TSV="${SEX_TSV:-${REPO_DIR}/panel_sex.tsv}"
EXCLUDE="${EXCLUDE-${REPO_DIR}/qc_excluded_birds.tsv}"
EXF="${EXCLUDE:-${OUT}/.no_exclusions}"     # awk needs a non-empty file even when nothing is excluded
LEPC_Z="${LEPC_Z:-NW_026294758.1,NW_026294813.1}"
TOOLS_ENV="${TOOLS_ENV:-${CLUSTER_SCRATCH}/GROUSE/z_markers/env}"   # samtools, python with numpy/pandas/matplotlib
HMM_ENV="${HMM_ENV:-${OUT}/env}"                                    # ancestry_hmm (bioconda)
DAF="${DAF:-0.6}"; THIN="${THIN:-10000}"   # more, evenly spread sites: at dAF >= 0.8 autosomal sites were sparse and clustered
RATE="${RATE:-3.11}"           # cM/Mb, uniform: chicken sex-averaged genome average, 3,176 cM over 1,020 Mb (Groenen et al. 2009)
RATE_TABLE="${RATE_TABLE:-}"   # optional scaffold<TAB>cM/Mb (e.g. chicken rate of each scaffold's homologous chromosome)
NE="${NE:-10000}"; ERR="${ERR:-0.01}"; TIME0="${TIME0:-20}"; BOOT="${BOOT:-0}"
GROUPS_RUN="${GROUPS_RUN:-control:LEPC_ctrl edge:LEPC_edge contact:LEPC_sym}"
QOS="${QOS:-normal}"; ACCOUNT="${ACCOUNT:-dewoody}"; PARTITION="${PARTITION:-cpu}"
BIRDS_PER_TASK="${BIRDS_PER_TASK:-8}"

die() { echo "ERROR: $*" >&2; exit 1; }
log() { echo "[$(date '+%F %T')] $*"; }
tools() { export PATH="${TOOLS_ENV}/bin:${PATH}"; }
REF() { cat "${HI}/ref.path"; }

stage_setup() {
    rm -rf "${OUT}/counts" "${OUT}/runs"          # new sites invalidate earlier read counts and runs
    mkdir -p "$OUT"/{counts,runs,logs}; echo "#" > "${OUT}/.no_exclusions"
    if [[ ! -x "${HMM_ENV}/bin/ancestry_hmm" ]]; then
        module load conda >/dev/null 2>&1 || module load anaconda >/dev/null 2>&1 || true
        CONDA_PKGS_DIRS="${OUT}/pkgs" "$(command -v mamba || command -v conda)" create -y -p "$HMM_ENV" \
            --override-channels -c conda-forge -c bioconda ancestry_hmm
    fi
    "${HMM_ENV}/bin/ancestry_hmm" 2>&1 | head -3 || true
    tools
    [[ -d "${HI}/red" && -d "${HI}/gl" ]] || die "no hybrid_index output in $HI"
    # hold out CONTROL_N Texas/New Mexico LEPC from the reference (fixed seed) as a negative control
    python3 - "$PANEL_GROUPS" "$GROUPS_TSV" "$CONTROL_N" "$EXF" <<'PYC'
import sys, random
src, dst, n, excl = sys.argv[1], sys.argv[2], int(sys.argv[3]), sys.argv[4]
ex = {l.split("\t")[0] for l in open(excl) if l.strip() and not l.startswith("#")}
rows = [l.rstrip("\n").split("\t") for l in open(src)]
h = rows[0]; ig, ist = h.index("grp"), h.index("STATE_NAME")
pool = [r[0] for r in rows[1:] if r[ig] == "LEPC_allo" and r[ist] in ("Texas", "New Mexico") and r[0] not in ex]
ctrl = set(random.Random(1).sample(sorted(pool), min(n, len(pool))))
with open(dst, "w") as o:
    for r in rows:
        if r[0] in ctrl:
            r[ig] = "LEPC_ctrl"
        o.write("\t".join(r) + "\n")
print(f"[setup] {len(ctrl)} Texas/New Mexico LEPC held out as LEPC_ctrl")
PYC
    python3 "$PY" sites --red "${HI}/red" --mafs "${HI}/gl" --groups "$GROUPS_TSV" ${EXCLUDE:+--exclude "$EXCLUDE"} \
        --z "$LEPC_Z" --daf "$DAF" --thin "$THIN" --out "${OUT}/sites"
    # birds to count: every admixed group, QC-passed
    local pat; pat=$(for g in $GROUPS_RUN; do echo "${g#*:}"; done | paste -sd'|')
    awk -F'\t' -v p="^(${pat})$" 'NR==FNR{if ($1 !~ /^#/) x[$1]=1; next} FNR>1 && $4 ~ p && !($1 in x){print $1"\t"$NF}' \
        "$EXF" "$GROUPS_TSV" > "${OUT}/birds.tsv"
    log "$(wc -l < "${OUT}/birds.tsv") admixed birds to count; sites in ${OUT}/sites.tsv"
}

stage_count() {
    tools
    local i=${SLURM_ARRAY_TASK_ID:?} ref; ref=$(REF)
    awk -v n="$BIRDS_PER_TASK" -v i="$i" 'NR > i*n && NR <= (i+1)*n' "${OUT}/birds.tsv" | while IFS=$'\t' read -r id cram; do
        [[ -s "${OUT}/counts/${id}.counts" ]] && continue
        samtools mpileup -l "${OUT}/sites.bed" -f "$ref" -q 20 -Q 20 -B "$cram" 2>/dev/null \
            | python3 "$PY" count --sites "${OUT}/sites.tsv" --out "${OUT}/counts/${id}.counts.tmp" \
            && mv "${OUT}/counts/${id}.counts.tmp" "${OUT}/counts/${id}.counts"
    done
}

stage_run() {
    tools
    local runs=() g c
    for g in $GROUPS_RUN; do for c in A Z; do runs+=("${g}:${c}"); done; done
    local r=${runs[${SLURM_ARRAY_TASK_ID:?}]}
    local name=${r%%:*} rest=${r#*:}; local grp=${rest%%:*} cls=${rest##*:}
    local tag="$cls" tspec="-${TIME0}"
    if [[ "${ZFIX:-0}" == 1 ]]; then
        [[ "$cls" == Z ]] || die "ZFIX runs are Z only"
        tspec=$(sed -n '/^optimum/,$p' "${OUT}/runs/${name}_A/ancestry_hmm.log" | awk '$1==1{print $2; exit}')
        [[ -n "$tspec" ]] || die "no autosomal pulse time in ${OUT}/runs/${name}_A/ancestry_hmm.log"
        tag=Zfix; log "$name Z: pulse time fixed at the autosomal estimate, $tspec generations"
    fi
    local dir="${OUT}/runs/${name}_${tag}"; mkdir -p "$dir"
    python3 "$PY" build --sites "${OUT}/sites.tsv" --counts "${OUT}/counts" --groups "$GROUPS_TSV" --sex "$SEX_TSV" \
        ${EXCLUDE:+--exclude "$EXCLUDE"} --grp "$grp" --cls "$cls" --rate-cm-mb "$RATE" \
        ${RATE_TABLE:+--rate-table "$RATE_TABLE"} --out "${dir}/data"
    # GRPC proportion = mean autosomal hybrid index of the group's (QC-passed) birds
    local m
    m=$(awk -F'\t' -v g="$grp" 'FNR==1{f++} f==1{if ($1 !~ /^#/) x[$1]=1; next}
        f==2{if (FNR==1){for(i=1;i<=NF;i++) gh[$i]=i} else if ($gh["grp"]==g) want[$1]=1; next}
        FNR==1{for(i=1;i<=NF;i++) h[$i]=i; next}
        ($h["ID"] in want) && $h["class"]=="autosomes" && !($h["ID"] in x){s+=$h["S"]; n++}
        END{v=(n ? s/n : 0); if (v < 0.001) v = 0.001; printf "%.4f", v}' \
        "$EXF" "$GROUPS_TSV" "${HI}/hybrid_index.tsv")
    local l; l=$(awk -v m="$m" 'BEGIN{printf "%.4f", 1-m}')
    log "$name $cls: GRPC proportion $m"
    cd "$dir"
    local bopt=""; [[ "$BOOT" -gt 0 ]] && bopt="-b $BOOT 1000"
    "${HMM_ENV}/bin/ancestry_hmm" -i data.input -s data.samples -a 2 "$l" "$m" \
        -p 0 100000 "$l" -p 1 "$tspec" "$m" --ne "$NE" -e "$ERR" $bopt \
        > ancestry_hmm.log 2>&1 || { tail -20 ancestry_hmm.log; die "ancestry_hmm failed: $name $cls"; }
    sed -n '/optimal model found/,$p' ancestry_hmm.log | head -8 || true     # estimated pulse time (generations)
}

stage_summary() {
    tools
    python3 "$PY" summarize --runs "${OUT}"/runs/*/ --z "$LEPC_Z" --out "${OUT}/ancestry_hmm" | tee "${OUT}/ancestry_hmm_summary.txt"
}

stage_submit() {
    [[ -s "${OUT}/birds.tsv" ]] || die "run: bash $0 setup"
    local nb nt nr J1 J2 J3
    nb=$(wc -l < "${OUT}/birds.tsv"); nt=$(( (nb + BIRDS_PER_TASK - 1) / BIRDS_PER_TASK ))
    nr=$(( $(wc -w <<< "$GROUPS_RUN") * 2 ))
    J1=$(sbatch --parsable -A "$ACCOUNT" -p "$PARTITION" -q "$QOS" -t 08:00:00 --array=0-$((nt - 1)) -c 2 --mem=8G \
        --export=ALL -J ahmm_count -o "${OUT}/logs/%x_%A_%a.out" "$SCRIPT" count) || die "count refused"
    J2=$(sbatch --parsable -A "$ACCOUNT" -p "$PARTITION" -q "$QOS" -t 1-00:00:00 --array=0-$((nr - 1)) -c 4 --mem=32G \
        --dependency=afterok:"$J1" --export=ALL -J ahmm_run -o "${OUT}/logs/%x_%A_%a.out" "$SCRIPT" run) || die "run refused"
    J3=$(sbatch --parsable -A "$ACCOUNT" -p "$PARTITION" -q "$QOS" -t 02:00:00 -c 2 --mem=32G \
        --dependency=afterok:"$J2" --export=ALL -J ahmm_summary -o "${OUT}/logs/%x_%j.out" "$SCRIPT" summary) || die "summary refused"
    log "counts $J1 ($nt tasks, $nb birds) -> runs $J2 ($nr) -> summary $J3"
}

stage_submit_zfix() {
    local nr J2 J3 idx
    nr=$(wc -w <<< "$GROUPS_RUN"); idx=$(seq -s, 1 2 $((2 * nr - 1)))      # the Z task of each group
    J2=$(sbatch --parsable -A "$ACCOUNT" -p "$PARTITION" -q "$QOS" -t 1-00:00:00 --array="$idx" -c 4 --mem=32G \
        --export=ALL,ZFIX=1 -J ahmm_zfix -o "${OUT}/logs/%x_%A_%a.out" "$SCRIPT" run) || die "run refused"
    J3=$(sbatch --parsable -A "$ACCOUNT" -p "$PARTITION" -q "$QOS" -t 02:00:00 -c 2 --mem=32G \
        --dependency=afterok:"$J2" --export=ALL -J ahmm_summary -o "${OUT}/logs/%x_%j.out" "$SCRIPT" summary) || die "summary refused"
    log "Z runs with fixed pulse time $J2 (tasks $idx) -> summary $J3"
}

stage_status() {
    echo "birds counted: $(ls "${OUT}"/counts/*.counts 2>/dev/null | wc -l) / $(wc -l < "${OUT}/birds.tsv")"
    for d in "${OUT}"/runs/*/; do [[ -d "$d" ]] && echo "$(basename "$d"): $(ls "$d"/*.posterior 2>/dev/null | wc -l) posterior files"; done
    squeue -u "$USER" -n ahmm_count,ahmm_run,ahmm_zfix,ahmm_summary 2>/dev/null | head
}

case "${1:-}" in
    setup) stage_setup ;; count) stage_count ;; run) stage_run ;; summary) stage_summary ;;
    submit) stage_submit ;; submit_zfix) stage_submit_zfix ;; status) stage_status ;;
    *) sed -n '2,25p' "$0" | sed 's/^# \{0,1\}//'; exit 1 ;;
esac
