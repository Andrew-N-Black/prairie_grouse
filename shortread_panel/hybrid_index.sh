#!/bin/bash
# =============================================================================
# HYBRID INDEX AND INTERSPECIFIC HETEROZYGOSITY, AUTOSOMES vs Z (Gautschi)
# Unnumbered (shortread_panel/). Uses the 506-bird low-coverage panel (crams_4.66x, aligned to the
# LEPC reference GCF_026119805.1) listed in panel_groups.tsv, with sex from panel_sex.tsv.
#
# By default the 2026 birds in panel_exclude_2026.txt are left out (held for the hybrid-zone paper):
# only the 24 allopatric 2026 GRPC are used, as parental reference birds (the 433-bird panel has 4).
# The 2026 STGR are among those held back, so only LEPC-GRPC is scored.
#
#   1. gl        (array over 5-Mb chunks of every autosomal scaffold >= MIN_SCAF and of the Z scaffolds)
#                cut each bird's reads for the chunk to node-local BAMs (ANGSD fails to decode some CRAM
#                slices; the Z uses step 22's Z-only BAMs when present), ANGSD genotype likelihoods for
#                all birds, then hybrid_index.py reduce keeps sites where an allopatric reference pair
#                differs in allele frequency by >= PRE_DAF (GLs of every bird kept, beagle deleted).
#   2. estimate  hybrid_index.py estimate: markers with dAF >= DAF between allopatric reference groups,
#                one per THIN bp; per bird and pair, ML hybrid index S and interspecific heterozygosity
#                H on the autosomes and the Z (females: S only, one Z), reference birds leave-one-out,
#                1-Mb block bootstrap; then triangle plots.
#
#   bash hybrid_index.sh setup      # login node: reference, chunk list, bam/cram list, python packages
#   bash hybrid_index.sh submit     # gl array -> estimate (+ plot)
#   bash hybrid_index.sh status
#   bash hybrid_index.sh submit-estimate   # rerun only step 2 (e.g. with DAF=0.9 or THIN=10000)
#
# Output ($OUT): hybrid_index.tsv (one row per bird x pair x chromosome class), hybrid_index.markers.txt,
#   hybrid_index_<pair>.png/pdf
# =============================================================================
#SBATCH -A dewoody
#SBATCH -p cpu
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=blackan@purdue.edu
set -euo pipefail
unset LD_PRELOAD

REPO_DIR="${SLURM_SUBMIT_DIR:-$PWD}"                    # run from shortread_panel/
SCRIPT="${REPO_DIR}/hybrid_index.sh"
PY="${REPO_DIR}/hybrid_index.py"
PANEL_TSV="${PANEL_TSV:-${REPO_DIR}/panel_groups.tsv}"
EXCLUDE="${EXCLUDE-${REPO_DIR}/panel_exclude_2026.txt}"   # birds held back (first column = ID); EXCLUDE= uses all
SEX_TSV="${SEX_TSV:-${REPO_DIR}/panel_sex.tsv}"
OUT="${OUT:-${CLUSTER_SCRATCH}/GROUSE/hybrid_index}"
GROUPS_TSV="${OUT}/groups.tsv"                           # panel minus EXCLUDE, written by setup
ENV_PREFIX="${ENV_PREFIX:-${CLUSTER_SCRATCH}/GROUSE/z_markers/env}"   # step 22's env (angsd, samtools)
ZBAM_DIR="${ZBAM_DIR:-${CLUSTER_SCRATCH}/GROUSE/z_markers/zbam}"      # step 22's Z-only BAMs
REF="${REF:-}"                                           # auto-detected from a CRAM header
LEPC_Z="${LEPC_Z:-NW_026294758.1,NW_026294813.1}"       # Z scaffolds of the LEPC reference (step 22)
MIN_SCAF="${MIN_SCAF:-1000000}"     # autosomal scaffolds shorter than this are skipped
CHUNK="${CHUNK:-5000000}"
PRE_DAF="${PRE_DAF:-0.5}"           # reduce: keep sites where any allopatric pair differs by this much
DAF="${DAF:-0.8}"                   # estimate: marker threshold for the pair being scored
THIN="${THIN:-5000}"                # one marker per THIN bp
BOOT="${BOOT:-100}"
PAIRS="${PAIRS:-LEPC-GRPC}"       # GRPC-STGR needs the 2026 STGR (EXCLUDE= PAIRS=LEPC-GRPC,GRPC-STGR)
MIN_IND_FRAC="${MIN_IND_FRAC:-0.5}"

ACCOUNT="${ACCOUNT:-dewoody}"; PARTITION="${PARTITION:-cpu}"
QOS="${QOS:-standby}"; WALL="${WALL:-04:00:00}"          # standby: idle nodes, 4 h limit
EST_QOS="${EST_QOS:-normal}"; EST_WALL="${EST_WALL:-1-00:00:00}"
MAXPAR="${MAXPAR:-40}"                                   # concurrent gl tasks
THREADS="${SLURM_CPUS_PER_TASK:-16}"

die() { echo "ERROR: $*" >&2; exit 1; }
log() { echo "[$(date '+%F %T')] $*"; }
use_env() { [[ -x "${ENV_PREFIX}/bin/angsd" ]] || die "no angsd in ${ENV_PREFIX} (run 22_z_markers.sh setup)"; export PATH="${ENV_PREFIX}/bin:${PATH}"; }
find_ref() {
    [[ -n "$REF" ]] && return
    [[ -s "${OUT}/ref.path" ]] && { REF="$(cat "${OUT}/ref.path")"; return; }
    local cram; cram="$(awk -F'\t' 'NR==2{print $NF}' "$PANEL_TSV")"
    REF="$(samtools view -H "$cram" | grep -m1 -o 'UR:[^[:space:]]*' | sed 's/^UR://; s#^file://##')"
    [[ -s "$REF" ]] || die "reference not found from CRAM header; set REF="
}

stage_setup() {
    use_env; mkdir -p "$OUT"/{gl,red,logs}
    find_ref; echo "$REF" > "${OUT}/ref.path"; log "reference: $REF"
    [[ -s "${REF}.fai" ]] || samtools faidx "$REF"
    python3 -c "import numpy, pandas, matplotlib" 2>/dev/null || {
        log "adding numpy/pandas/matplotlib to ${ENV_PREFIX}"
        module load conda >/dev/null 2>&1 || module load anaconda >/dev/null 2>&1 || true
        CONDA_PKGS_DIRS="${CLUSTER_SCRATCH}/GROUSE/z_markers/pkgs" \
            "$(command -v mamba || command -v conda)" install -y -p "$ENV_PREFIX" -c conda-forge numpy pandas matplotlib; }
    # panel minus the held-back birds; bird order = that file's order; one CRAM per bird
    if [[ -n "$EXCLUDE" ]]; then
        [[ -s "$EXCLUDE" ]] || die "EXCLUDE=$EXCLUDE not found"
        awk -F'\t' 'NR==FNR{if ($1 !~ /^#/) x[$1]=1; next} FNR==1 || !($1 in x)' "$EXCLUDE" "$PANEL_TSV" > "$GROUPS_TSV"
    else
        cp "$PANEL_TSV" "$GROUPS_TSV"
    fi
    log "birds used: $(($(wc -l < "$GROUPS_TSV") - 1)) of $(($(wc -l < "$PANEL_TSV") - 1)) ($(awk -F'\t' 'NR>1{print $4}' "$GROUPS_TSV" | sort | uniq -c | tr -s ' ' | tr '\n' ';'))"
    REFS=$(awk -F'\t' 'NR>1 && $4 ~ /_allo$/ {print $4}' "$GROUPS_TSV" | sort -u | paste -sd,); echo "$REFS" > "${OUT}/refs.txt"
    awk -F'\t' 'NR>1{print $1}' "$GROUPS_TSV" > "${OUT}/ids.txt"
    awk -F'\t' 'NR>1{print $1"\t"$NF}' "$GROUPS_TSV" > "${OUT}/crams.tsv"
    while IFS=$'\t' read -r id c; do [[ -s "$c" ]] || die "missing CRAM for $id: $c"; done < "${OUT}/crams.tsv"
    # chunks: autosomal scaffolds >= MIN_SCAF, plus the Z scaffolds (any length)
    awk -v z="$LEPC_Z" -v m="$MIN_SCAF" -v c="$CHUNK" 'BEGIN{n=split(z,a,","); for(i=1;i<=n;i++) Z[a[i]]=1}
        { if (!($1 in Z) && $2 < m) next
          for (s=1; s<=$2; s+=c) { e=s+c-1; if (e>$2) e=$2; print $1":"s"-"e"\t"(($1 in Z)?"Z":"A") } }' \
        "${REF}.fai" > "${OUT}/chunks.tsv"
    log "$(wc -l < "${OUT}/ids.txt") birds; $(cut -f2 "${OUT}/chunks.tsv" | sort | uniq -c | tr '\n' ' ')chunks (A = autosomal, Z)"
    log "autosomal sequence analysed: $(awk -v z="$LEPC_Z" -v m="$MIN_SCAF" 'BEGIN{n=split(z,a,","); for(i=1;i<=n;i++) Z[a[i]]=1} !($1 in Z) && $2>=m {s+=$2} END{printf "%.0f Mb", s/1e6}' "${REF}.fai")"
}

stage_gl() {
    use_env; find_ref
    local line region cls tag tmp nb
    line="$(sed -n "$((SLURM_ARRAY_TASK_ID + 1))p" "${OUT}/chunks.tsv")"
    region="${line%%$'\t'*}"; cls="${line##*$'\t'}"; tag="${region//[:]/_}"
    [[ -s "${OUT}/red/${tag}.npz" ]] && { log "$tag done"; return 0; }
    tmp="${OUT}/tmp/${tag}"; rm -rf "$tmp"; mkdir -p "$tmp"     # ~5 GB of chunk BAMs; scratch, not node /tmp
    trap 'rm -rf "$tmp"' EXIT
    : > "${tmp}/bams.txt"
    export REF tmp region
    one() {   # $1 = id, $2 = cram
        local b="${tmp}/$1.bam"
        samtools view -b -T "$REF" -o "$b" "$2" "$region" && samtools index "$b" || { echo "FAILED $1" >&2; return 1; }
    }
    export -f one
    if [[ "$cls" == "Z" && -d "$ZBAM_DIR" ]]; then
        while IFS=$'\t' read -r id c; do
            b="${ZBAM_DIR}/$(basename "$c" .cram).z.bam"
            [[ -s "${b}.bai" ]] || die "no Z-only BAM for $id ($b); unset ZBAM_DIR to cut from CRAMs"
            echo "$b" >> "${tmp}/bams.txt"
        done < "${OUT}/crams.tsv"
    else
        xargs -a "${OUT}/crams.tsv" -P "$THREADS" -L1 bash -c 'one "$0" "$1"' || die "BAM cutting failed for $region"
        while IFS=$'\t' read -r id c; do echo "${tmp}/${id}.bam" >> "${tmp}/bams.txt"; done < "${OUT}/crams.tsv"
    fi
    nb=$(wc -l < "${tmp}/bams.txt")
    angsd -bam "${tmp}/bams.txt" -ref "$REF" -r "$region" -out "${tmp}/${tag}" -P 4 \
        -uniqueOnly 1 -remove_bads 1 -only_proper_pairs 1 -minMapQ 30 -minQ 20 -C 50 -baq 1 \
        -GL 2 -doGlf 2 -doMajorMinor 1 -doMaf 1 -SNP_pval 1e-6 -minMaf 0.02 \
        -doCounts 1 -minInd "$(awk -v n="$nb" -v f="$MIN_IND_FRAC" 'BEGIN{printf "%d", n*f}')" -setMaxDepth $((nb * 15)) \
        > "${tmp}/${tag}.angsd.log" 2>&1 || { tail -20 "${tmp}/${tag}.angsd.log"; die "angsd failed: $region"; }
    python3 "$PY" reduce --beagle "${tmp}/${tag}.beagle.gz" --ids "${OUT}/ids.txt" --groups "$GROUPS_TSV" \
        --refs "$(cat "${OUT}/refs.txt")" --pre-daf "$PRE_DAF" --out "${tmp}/${tag}.npz"
    cp "${tmp}/${tag}.mafs.gz" "${OUT}/gl/" && mv "${tmp}/${tag}.npz" "${OUT}/red/${tag}.npz"
    log "$region ($cls): done"
}

stage_estimate() {
    use_env
    python3 "$PY" estimate --indir "${OUT}/red" --groups "$GROUPS_TSV" --sex "$SEX_TSV" --z "$LEPC_Z" \
        --pairs "$PAIRS" --daf "$DAF" --thin "$THIN" --boot "$BOOT" --threads "$THREADS" --out "${OUT}/hybrid_index.tsv"
    python3 "$PY" plot --tsv "${OUT}/hybrid_index.tsv" --prefix "${OUT}/hybrid_index"
    log "wrote ${OUT}/hybrid_index.tsv and plots"
}

submit_estimate() {
    local dep="${1:-}"
    sbatch --parsable -A "$ACCOUNT" -p "$PARTITION" -q "$EST_QOS" -t "$EST_WALL" --export=ALL ${dep:+--dependency=afterok:$dep} \
        -c 64 --mem=200G -J hi_estimate -o "${OUT}/logs/%x_%j.out" "$SCRIPT" estimate
}

stage_submit() {
    [[ -s "${OUT}/chunks.tsv" ]] || die "run: bash $0 setup"
    local n J1 J2; n=$(wc -l < "${OUT}/chunks.tsv")
    J1=$(sbatch --parsable -A "$ACCOUNT" -p "$PARTITION" -q "$QOS" -t "$WALL" --export=ALL \
        --array=0-$((n - 1))%"$MAXPAR" -c 16 --mem=32G -J hi_gl -o "${OUT}/logs/%x_%A_%a.out" "$SCRIPT" gl) || die "gl refused"
    J2=$(submit_estimate "$J1") || die "estimate refused"
    log "gl array $J1 ($n chunks) -> estimate $J2"
    log "standby tasks that hit the 4-h limit can be resubmitted with the same command (finished chunks are skipped);"
    log "the estimate job then needs: bash $0 submit-estimate"
}

stage_status() {
    local n d; n=$(wc -l < "${OUT}/chunks.tsv"); d=$(ls "${OUT}"/red/*.npz 2>/dev/null | wc -l)
    echo "chunks done: $d / $n"
    squeue -u "$USER" -n hi_gl,hi_estimate -o "%.12i %.12j %.8T %.10M %R" 2>/dev/null | head -20
    [[ -s "${OUT}/hybrid_index.markers.txt" ]] && cat "${OUT}/hybrid_index.markers.txt"
    grep -l -E 'ERROR|FAILED|CANCELLED|TIME LIMIT' "${OUT}"/logs/hi_gl_*.out 2>/dev/null | head || true
}

case "${1:-}" in
    setup) stage_setup ;; gl) stage_gl ;; estimate) stage_estimate ;;
    submit) stage_submit ;; submit-estimate) J=$(submit_estimate) && log "estimate $J" ;; status) stage_status ;;
    *) sed -n '2,28p' "$0" | sed 's/^# \{0,1\}//'; exit 1 ;;
esac
