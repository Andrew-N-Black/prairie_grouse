#!/bin/bash
# =============================================================================
# STEP 14: SPECIES-DIAGNOSTIC Z-CHROMOSOME MARKERS (Gautschi)
#
# Two-stage design:
#   1. DISCOVER in the 46 haplotype assemblies (step 05/06 output): align each haplotype's chr_Z
#      to the LEPC reference Z scaffolds (the coordinates of the 506-bird panel), call SNPs and
#      indels per haplotype (minimap2 + paftools call), and keep sites where every called
#      haplotype of one species carries one allele and every called haplotype of another species
#      a different one. Then check each candidate's ~300-bp window: single copy in the LEPC
#      genome, no match to any female W (gametologs would also amplify from the W), few repeats,
#      and design primers around it (Primer3; other variable sites are excluded from priming).
#   2. CONFIRM in the 506-bird low-coverage panel: allele frequencies at the candidate SNPs in
#      each species x context group (ANGSD genotype likelihoods, -sites), so "fixed" means fixed
#      across hundreds of birds, not 7-19 assembled Z copies. Sympatric groups show how often
#      the other species' allele turns up in contact zones.
#
#   bash 14_z_markers.sh setup       # once, login node: conda env, reference Z, W, bam lists
#   bash 14_z_markers.sh submit      # align (array) -> discover -> panel (array) -> summarize
#                                    # all on dewoody standby (QOS=standby, 4 h); QOS=normal WALL=1-00:00:00 to change
#   bash 14_z_markers.sh status
#   bash 14_z_markers.sh select        # re-pick the panel SNPs (MIN_CLEAN, BIN, MAX_MASKED) on the login node
#   bash 14_z_markers.sh submit-panel  # rerun only the ANGSD panel check + summary
#   bash 14_z_markers.sh submit-zbam   # cut Z-only BAMs from the CRAMs first (ANGSD CRAM decode errors), then panel + summary
#
# Output: ${OUT_DIR}/z_markers.final.tsv, one row per candidate: alleles per species in the
#   assemblies, haplotypes called, flank cleanliness, genome copies, W hits, repeat fraction,
#   primer pair, ANGSD alt-allele frequency and nInd per group (LEPC/GRPC/STGR x allo/sym, HYB),
#   allele-frequency difference between allopatric groups for each pair, and
#   panel_confirmed_pairs (difference >= MIN_DAF). Sorted: most pairs confirmed first, indels first.
#
# NOTES
#   * Indels (best for gel/field assays) cannot be checked by ANGSD; they carry the assembly
#     evidence only and need the lab validation set. A diagnostic SNP inside the same window
#     can stand in for them in the panel table.
#   * Females are hemizygous for the Z; ANGSD treats them as diploid homozygotes. That does not
#     affect whether an allele is fixed, only slightly how frequencies are weighted.
#   * A Z marker reads the bird's paternal Z: an F1 female looks like her father's species.
#     Pair the panel with CHD1 sexing, mtDNA and a few autosomal markers for hybrid screening.
# =============================================================================
#SBATCH -A dewoody
#SBATCH -p cpu
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=blackan@purdue.edu

set -euo pipefail
unset LD_PRELOAD

REPO_DIR="${SLURM_SUBMIT_DIR:-$PWD}"
SCRIPT="${REPO_DIR}/14_z_markers.sh"
PYHELP="${REPO_DIR}/14_z_markers.py"
GROUPS_TSV="${GROUPS_TSV:-${REPO_DIR}/panel_groups.tsv}"   # ID, species, context, group, cram

ASM="${ASM:-${CLUSTER_SCRATCH}/GROUSE/grouse_asm}"
ASM_GLOB="${ASM_GLOB:-${ASM}/final_ptarmigan/*_hap[12].pseudo_chr.fasta}"
ASM_Z="${ASM_Z:-chr_Z}"
OUT_DIR="${OUT_DIR:-${CLUSTER_SCRATCH}/GROUSE/z_markers}"
ENV_PREFIX="${ENV_PREFIX:-${OUT_DIR}/env}"
LEPC_REF="${LEPC_REF:-}"                       # FASTA the panel CRAMs were aligned to (auto-detected)
LEPC_Z="${LEPC_Z:-NW_026294758.1,NW_026294813.1}"   # Z scaffolds of the LEPC reference

ACCOUNT="${ACCOUNT:-dewoody}"
PARTITION="${PARTITION:-cpu}"
QOS="${QOS:-standby}"                    # dewoody normal queue is busy; standby = idle nodes, 4 h limit
WALL="${WALL:-04:00:00}"
THREADS="${SLURM_CPUS_PER_TASK:-16}"

MIN_CALLED="${MIN_CALLED:-0.8}"     # fraction of a species' haplotypes that must be called
MIN_COPIES="${MIN_COPIES:-4}"       # and at least this many Z copies per species
FLANK="${FLANK:-150}"               # window = variant +/- FLANK bp (primers are placed in it)
MAX_MASKED="${MAX_MASKED:-0.2}"     # max soft-masked (repeat) fraction of the window
MIN_CLEAN="${MIN_CLEAN:-20}"        # no other variable site within this many bp of the marker
BIN="${BIN:-20000}"                 # panel: keep the cleanest marker per BIN bp per species pair
MIN_DAF="${MIN_DAF:-0.95}"          # allele-frequency difference between allopatric groups
MIN_IND_FRAC="${MIN_IND_FRAC:-0.5}" # ANGSD -minInd as a fraction of the group
# PANEL_RF: regions ANGSD reads. Default = one region per candidate SNP (fast). If a group fails
#   with "Failure to decode slice" (CRAM random access), use the whole Z instead:
#   PANEL_RF=${OUT_DIR}/ref/z_regions.txt  (slower; output is still limited to the candidate SNPs)
PRODUCT="${PRODUCT:-80-200}"        # short products for feathers / degraded field samples

CALL_DIR="${OUT_DIR}/calls"; DISC="${OUT_DIR}/discovery"; PANEL="${OUT_DIR}/panel"
MODE="${1:-}"
die() { echo "ERROR: $*" >&2; exit 1; }
log() { echo "[$(date '+%F %T')] $*"; }
use_env() { [[ -x "${ENV_PREFIX}/bin/angsd" ]] || die "run: bash $0 setup"; export PATH="${ENV_PREFIX}/bin:${PATH}"; }

find_ref() {
    [[ -n "$LEPC_REF" ]] && { [[ -s "$LEPC_REF" ]] || die "LEPC_REF=$LEPC_REF not found"; return; }
    [[ -s "${OUT_DIR}/ref/lepc_ref.path" ]] && { LEPC_REF="$(cat "${OUT_DIR}/ref/lepc_ref.path")"; return; }
    local cram ur
    cram="$(awk -F'\t' 'NR==2{print $NF}' "$GROUPS_TSV")"
    ur="$(samtools view -H "$cram" 2>/dev/null | grep -m1 -o 'UR:[^[:space:]]*' | sed 's/^UR://; s#^file://##' || true)"
    [[ -s "$ur" ]] || die "could not find the panel reference from the CRAM header (UR: '$ur'); set LEPC_REF=/path/to/LEPC.fna"
    LEPC_REF="$ur"
}

# ----------------------------------------------------------------------------- setup
stage_setup() {
    mkdir -p "$OUT_DIR"/{ref,calls,discovery,panel/bamlists} "${REPO_DIR}/logs"
    [[ -s "$GROUPS_TSV" ]] || die "missing $GROUPS_TSV"
    if [[ ! -x "${ENV_PREFIX}/bin/angsd" ]]; then
        module load conda >/dev/null 2>&1 || module load anaconda >/dev/null 2>&1 || true
        command -v conda >/dev/null || die "no conda"
        export CONDA_PKGS_DIRS="${OUT_DIR}/pkgs"
        "$(command -v mamba || command -v conda)" create -y -p "$ENV_PREFIX" --override-channels \
            -c conda-forge -c bioconda minimap2 samtools htslib angsd primer3 python=3.10
    fi
    use_env
    command -v paftools.js >/dev/null || die "paftools.js missing from the env (bioconda minimap2 ships it)"
    find_ref; echo "$LEPC_REF" > "${OUT_DIR}/ref/lepc_ref.path"; log "LEPC reference: $LEPC_REF"
    [[ -s "${LEPC_REF}.fai" ]] || samtools faidx "$LEPC_REF"
    samtools faidx "$LEPC_REF" ${LEPC_Z//,/ } > "${OUT_DIR}/ref/lepcZ.fa"
    tr ',' '\n' <<< "$LEPC_Z" > "${OUT_DIR}/ref/z_regions.txt"
    # haplotypes with a Z (W-bearing female haplotypes have almost none) and female W sequence
    : > "${OUT_DIR}/ref/haplotypes.txt"; : > "${OUT_DIR}/ref/female_W.fa"
    for f in $ASM_GLOB; do
        [[ -s "${f}.fai" ]] || samtools faidx "$f"
        local h zl wl; h="$(basename "$f" .pseudo_chr.fasta)"
        zl=$(awk -v z="$ASM_Z" '$1==z{print $2}' "${f}.fai"); wl=$(awk '$1=="chr_W"{print $2}' "${f}.fai")
        [[ "${zl:-0}" -ge 20000000 ]] && echo -e "${h}\t${f}" >> "${OUT_DIR}/ref/haplotypes.txt"
        [[ "${wl:-0}" -ge 1000000 ]] && samtools faidx "$f" chr_W | sed "s/^>.*/>${h}_W/" >> "${OUT_DIR}/ref/female_W.fa"
    done
    log "haplotypes with a Z: $(wc -l < "${OUT_DIR}/ref/haplotypes.txt"); female W sequences: $(grep -c '>' "${OUT_DIR}/ref/female_W.fa")"
    # one bam list per group (CRAMs)
    awk -F'\t' -v d="${PANEL}/bamlists" 'NR>1{print $NF > (d "/" $4 ".list")}' "$GROUPS_TSV"
    for l in "${PANEL}"/bamlists/*.list; do printf "  %-12s %d birds\n" "$(basename "$l" .list)" "$(wc -l < "$l")"; done
}

# ----------------------------------------------------------------------------- align (array)
stage_align() {
    use_env
    local t="${SLURM_ARRAY_TASK_ID:?array job}" h f
    read -r h f < <(sed -n "$((t + 1))p" "${OUT_DIR}/ref/haplotypes.txt")
    [[ -n "$h" ]] || { log "no haplotype for task $t"; return 0; }
    [[ -s "${CALL_DIR}/${h}.call.txt" ]] && { log "$h done"; return 0; }
    log "$h: chr_Z -> LEPC Z"
    samtools faidx "$f" "$ASM_Z" > "${CALL_DIR}/${h}.z.fa"
    minimap2 -c -x asm5 --cs -t "$THREADS" "${OUT_DIR}/ref/lepcZ.fa" "${CALL_DIR}/${h}.z.fa" \
        | sort -k6,6 -k8,8n > "${CALL_DIR}/${h}.paf"
    paftools.js call -l 10000 -L 50000 -q 5 "${CALL_DIR}/${h}.paf" > "${CALL_DIR}/${h}.call.tmp"
    mv "${CALL_DIR}/${h}.call.tmp" "${CALL_DIR}/${h}.call.txt"
    rm -f "${CALL_DIR}/${h}.z.fa"
    log "$h: $(grep -c '^V' "${CALL_DIR}/${h}.call.txt") variants; $(awk '$1=="R"{s+=$4-$3}END{printf "%.1f", s/1e6}' "${CALL_DIR}/${h}.call.txt") Mb called"
}

# ----------------------------------------------------------------------------- discover
stage_discover() {
    use_env; find_ref
    [[ -s "${OUT_DIR}/ref/lepc_genome.mmi" ]] || minimap2 -x sr -t "$THREADS" -d "${OUT_DIR}/ref/lepc_genome.mmi" "$LEPC_REF"
    log "merging $(ls "$CALL_DIR"/*.call.txt | wc -l) haplotypes"
    python3 "$PYHELP" discover --calls "$CALL_DIR" --out "$DISC" --min-called "$MIN_CALLED" \
        --min-copies "$MIN_COPIES" --flank "$FLANK"
    # window sequences from the soft-masked LEPC reference, named by candidate id
    awk '{printf "%s:%d-%d\t%s\n", $1, $2+1, $3, $4}' "${DISC}/windows.bed" > "${DISC}/windows.regions"
    cut -f1 "${DISC}/windows.regions" > "${DISC}/regions.txt"
    samtools faidx "$LEPC_REF" -r "${DISC}/regions.txt" \
        | awk 'NR==FNR{id[">"$1]=$2; next} /^>/{print ">" id[$1]; next} {print}' "${DISC}/windows.regions" - \
        > "${DISC}/amplicons.fa"
    log "copy number in the LEPC genome and matches to the female W"
    minimap2 -c -x sr -k 15 -w 5 --secondary=yes -N 50 -p 0.5 -t "$THREADS" "${OUT_DIR}/ref/lepc_genome.mmi" \
        "${DISC}/amplicons.fa" > "${DISC}/amplicons.genome.paf" 2>/dev/null
    if [[ -s "${OUT_DIR}/ref/female_W.fa" ]]; then
        minimap2 -c -x sr -k 15 -w 5 --secondary=yes -N 50 -p 0.5 -t "$THREADS" "${OUT_DIR}/ref/female_W.fa" \
            "${DISC}/amplicons.fa" > "${DISC}/amplicons.W.paf" 2>/dev/null
    else : > "${DISC}/amplicons.W.paf"; fi
    python3 "$PYHELP" hits --cand "${DISC}/candidates.raw.tsv" --amplicons "${DISC}/amplicons.fa" \
        --windows "${DISC}/windows.bed" --varpos "${DISC}/varpos.bed" --genome-paf "${DISC}/amplicons.genome.paf" \
        --w-paf "${DISC}/amplicons.W.paf" --primer3-in "${DISC}/primer3.in" --out "${DISC}/candidates.checked.tsv"
    cat > "${DISC}/primer3.settings" <<EOF
Primer3 File - http://primer3.org
P3_FILE_TYPE=settings

PRIMER_TASK=generic
PRIMER_PICK_LEFT_PRIMER=1
PRIMER_PICK_RIGHT_PRIMER=1
PRIMER_PRODUCT_SIZE_RANGE=${PRODUCT}
PRIMER_OPT_SIZE=20
PRIMER_MIN_SIZE=18
PRIMER_MAX_SIZE=25
PRIMER_OPT_TM=60.0
PRIMER_MIN_TM=57.0
PRIMER_MAX_TM=63.0
PRIMER_MIN_GC=35
PRIMER_MAX_GC=65
PRIMER_LOWERCASE_MASKING=1
PRIMER_NUM_RETURN=1
=
EOF
    primer3_core --p3_settings_file="${DISC}/primer3.settings" < "${DISC}/primer3.in" > "${DISC}/primer3.out" 2> "${DISC}/primer3.log" || true
    python3 "$PYHELP" primers --cand "${DISC}/candidates.checked.tsv" --p3out "${DISC}/primer3.out" \
        --out "${DISC}/candidates.primers.tsv"
    stage_select
    log "candidates: $(($(wc -l < "${DISC}/candidates.primers.tsv") - 1)); SNPs to the panel: $(wc -l < "${PANEL}/sites.txt")"
}

# ----------------------------------------------------------------------------- select (login node, seconds)
# SNPs for the panel: window checks passed, primer pair, clean flanks, then the cleanest one per
# BIN bp for each species pair. ANGSD needs fresh index files for both the sites and the reference.
stage_select() {
    use_env; find_ref
    python3 "$PYHELP" select --cand "${DISC}/candidates.primers.tsv" --out "${PANEL}/sites.txt" \
        --max-masked "$MAX_MASKED" --min-clean "$MIN_CLEAN" --bin "$BIN"
    rm -f "${PANEL}/sites.txt.bin" "${PANEL}/sites.txt.idx"
    angsd sites index "${PANEL}/sites.txt" > "${PANEL}/sites.index.log" 2>&1 \
        || { cat "${PANEL}/sites.index.log"; die "angsd sites index failed"; }
    awk '{print $1":"$2"-"$2}' "${PANEL}/sites.txt" > "${PANEL}/sites.regions"
    if [[ ! -s "${LEPC_REF}.fai" || ! "${LEPC_REF}.fai" -nt "$LEPC_REF" ]]; then
        log "re-indexing $LEPC_REF (ANGSD refuses a .fai that is not strictly newer than the FASTA)"
        samtools faidx "$LEPC_REF"
    fi
}

# ----------------------------------------------------------------------------- panel (array over groups)
stage_panel() {
    use_env; find_ref
    local ldir="${PANEL}/bamlists"; [[ -d "${PANEL}/bamlists_z" ]] && ldir="${PANEL}/bamlists_z"
    local lists=("${ldir}"/*.list) l g n
    l="${lists[${SLURM_ARRAY_TASK_ID:?array job}]}"; g="$(basename "$l" .list)"; n=$(wc -l < "$l")
    log "$g: $n birds (lists: $(basename "$ldir"))"
    [[ -s "${PANEL}/${g}.mafs.gz" ]] && { log "$g done"; return 0; }
    angsd -b "$l" -ref "$LEPC_REF" -sites "${PANEL}/sites.txt" -rf "${PANEL_RF:-${PANEL}/sites.regions}" \
        -GL 2 -doMajorMinor 3 -doMaf 1 -doCounts 1 \
        -minMapQ 20 -minQ 20 -remove_bads 1 -uniqueOnly 1 -only_proper_pairs 1 \
        -minInd "$(awk -v n="$n" -v f="$MIN_IND_FRAC" 'BEGIN{m=int(n*f); print (m<1?1:m)}')" \
        -nThreads "$THREADS" -out "${PANEL}/${g}.tmp" 2> "${PANEL}/${g}.angsd.log"
    [[ -s "${PANEL}/${g}.tmp.mafs.gz" ]] || { tail -20 "${PANEL}/${g}.angsd.log"; die "ANGSD wrote no output for $g"; }
    mv "${PANEL}/${g}.tmp.mafs.gz" "${PANEL}/${g}.mafs.gz"
    log "$g: $(($(zcat "${PANEL}/${g}.mafs.gz" | wc -l) - 1)) sites with data"
}

stage_summarize() {
    use_env
    python3 "$PYHELP" panel --cand "${DISC}/candidates.primers.tsv" --mafs-dir "$PANEL" \
        --out "${OUT_DIR}/z_markers.final.tsv" --min-daf "$MIN_DAF" \
        --sites "${PANEL}/sites.txt" --max-masked "$MAX_MASKED"
    log "final table: ${OUT_DIR}/z_markers.final.tsv"
}

stage_submit() {
    mkdir -p "${REPO_DIR}/logs"
    [[ -s "${OUT_DIR}/ref/haplotypes.txt" ]] || die "run: bash $0 setup"
    local NH NG J1 J2 J3 J4
    NH=$(wc -l < "${OUT_DIR}/ref/haplotypes.txt"); NG=$(ls "${PANEL}"/bamlists/*.list | wc -l)
    J1=$(sbatch --parsable -A "$ACCOUNT" -p "$PARTITION" -q "$QOS" -t "$WALL" --export=ALL --array=0-$((NH-1))%16 -c 16 --mem=32G \
        -J zm_align -o "${REPO_DIR}/logs/%x_%A_%a.out" "$SCRIPT" align) || die "align refused"
    J2=$(sbatch --parsable -A "$ACCOUNT" -p "$PARTITION" -q "$QOS" -t "$WALL" --export=ALL --dependency=afterok:"$J1" -c 16 --mem=64G \
        -J zm_discover -o "${REPO_DIR}/logs/%x_%j.out" "$SCRIPT" discover) || die "discover refused"
    J3=$(sbatch --parsable -A "$ACCOUNT" -p "$PARTITION" -q "$QOS" -t "$WALL" --export=ALL --dependency=afterok:"$J2" --array=0-$((NG-1)) \
        -c 16 --mem=32G -J zm_panel -o "${REPO_DIR}/logs/%x_%A_%a.out" "$SCRIPT" panel) || die "panel refused"
    J4=$(sbatch --parsable -A "$ACCOUNT" -p "$PARTITION" -q "$QOS" -t "$WALL" --export=ALL --dependency=afterok:"$J3" -c 2 --mem=8G \
        -J zm_summary -o "${REPO_DIR}/logs/%x_%j.out" "$SCRIPT" summarize) || die "summary refused"
    log "align $J1 ($NH haplotypes) -> discover $J2 -> panel $J3 ($NG groups) -> summarize $J4  [qos $QOS, $WALL]"
    echo "Standby jobs can be killed at the limit or preempted; every stage skips finished work, so resubmit the same stage:"
    echo "  sbatch -A $ACCOUNT -p $PARTITION -q $QOS -t $WALL --export=ALL --array=0-$((NG-1)) -c 16 --mem=32G -o logs/%x_%A_%a.out -J zm_panel $SCRIPT panel"
}

# ----------------------------------------------------------------------------- zbam (array over groups)
# ANGSD 0.940's CRAM reader fails ("Failure to decode slice", thousands of times) on the LEPC and
# GRPC CRAMs, silently dropping reads, although samtools reads them cleanly. So samtools cuts each
# bird's Z scaffolds into a small indexed BAM (~30 MB at 4.7x) and ANGSD reads those instead.
stage_zbam() {
    use_env; find_ref
    local lists=("${PANEL}"/bamlists/*.list) l g
    l="${lists[${SLURM_ARRAY_TASK_ID:?array job}]}"; g="$(basename "$l" .list)"
    mkdir -p "${OUT_DIR}/zbam" "${PANEL}/bamlists_z"
    export LEPC_REF OUT_DIR LEPC_Z
    one() {
        local c="$1" b; b="${OUT_DIR}/zbam/$(basename "$c" .cram).z.bam"
        [[ -s "${b}.bai" ]] && { echo "$b"; return 0; }
        samtools view -b -T "$LEPC_REF" -o "${b}.tmp" "$c" ${LEPC_Z//,/ } \
            && mv "${b}.tmp" "$b" && samtools index "$b" && echo "$b" \
            || { echo "FAILED $c" >&2; rm -f "${b}.tmp"; }
    }
    export -f one
    xargs -a "$l" -P "$THREADS" -I{} bash -c 'one "$1"' _ {} > "${PANEL}/bamlists_z/${g}.list.tmp"
    local n_in n_out; n_in=$(wc -l < "$l"); n_out=$(wc -l < "${PANEL}/bamlists_z/${g}.list.tmp")
    [[ "$n_in" -eq "$n_out" ]] || die "$g: only $n_out of $n_in Z BAMs made"
    sort "${PANEL}/bamlists_z/${g}.list.tmp" > "${PANEL}/bamlists_z/${g}.list" && rm "${PANEL}/bamlists_z/${g}.list.tmp"
    log "$g: $n_out Z-only BAMs"
}

stage_submit_zbam() {
    mkdir -p "${REPO_DIR}/logs"
    local NG J0 J3 J4
    NG=$(ls "${PANEL}"/bamlists/*.list | wc -l)
    J0=$(sbatch --parsable -A "$ACCOUNT" -p "$PARTITION" -q "$QOS" -t "$WALL" --export=ALL --array=0-$((NG-1)) \
        -c 16 --mem=16G -J zm_zbam -o "${REPO_DIR}/logs/%x_%A_%a.out" "$SCRIPT" zbam) || die "zbam refused"
    J3=$(sbatch --parsable -A "$ACCOUNT" -p "$PARTITION" -q "$QOS" -t "$WALL" --export=ALL --dependency=afterok:"$J0" \
        --array=0-$((NG-1)) -c 16 --mem=32G -J zm_panel -o "${REPO_DIR}/logs/%x_%A_%a.out" "$SCRIPT" panel) || die "panel refused"
    J4=$(sbatch --parsable -A "$ACCOUNT" -p "$PARTITION" -q "$QOS" -t "$WALL" --export=ALL --dependency=afterok:"$J3" -c 2 --mem=8G \
        -J zm_summary -o "${REPO_DIR}/logs/%x_%j.out" "$SCRIPT" summarize) || die "summary refused"
    log "zbam $J0 -> panel $J3 -> summarize $J4  [qos $QOS, $WALL]"
}

stage_submit_panel() {
    mkdir -p "${REPO_DIR}/logs"
    [[ -s "${PANEL}/sites.txt" ]] || die "run: bash $0 select"
    local NG J3 J4
    NG=$(ls "${PANEL}"/bamlists/*.list | wc -l)
    J3=$(sbatch --parsable -A "$ACCOUNT" -p "$PARTITION" -q "$QOS" -t "$WALL" --export=ALL --array=0-$((NG-1)) \
        -c 16 --mem=32G -J zm_panel -o "${REPO_DIR}/logs/%x_%A_%a.out" "$SCRIPT" panel) || die "panel refused"
    J4=$(sbatch --parsable -A "$ACCOUNT" -p "$PARTITION" -q "$QOS" -t "$WALL" --export=ALL --dependency=afterok:"$J3" -c 2 --mem=8G \
        -J zm_summary -o "${REPO_DIR}/logs/%x_%j.out" "$SCRIPT" summarize) || die "summary refused"
    log "panel $J3 ($NG groups) -> summarize $J4  [qos $QOS, $WALL]"
}

stage_status() {
    echo "calls:      $(ls "$CALL_DIR"/*.call.txt 2>/dev/null | wc -l) / $(wc -l < "${OUT_DIR}/ref/haplotypes.txt" 2>/dev/null || echo ?)"
    [[ -s "${DISC}/candidates.primers.tsv" ]] && echo "candidates: $(($(wc -l < "${DISC}/candidates.primers.tsv") - 1))"
    echo "panel:      $(ls "$PANEL"/*.mafs.gz 2>/dev/null | wc -l) groups done"
    [[ -s "${OUT_DIR}/z_markers.final.tsv" ]] && echo "final:      ${OUT_DIR}/z_markers.final.tsv"
    return 0
}

case "$MODE" in
    setup) stage_setup ;; align) stage_align ;; discover) stage_discover ;;
    select) stage_select ;; submit-panel) stage_submit_panel ;; zbam) stage_zbam ;; submit-zbam) stage_submit_zbam ;;
    panel) stage_panel ;; summarize) stage_summarize ;; submit) stage_submit ;; status) stage_status ;;
    *) sed -n '2,40p' "$0" | sed 's/^# \{0,1\}//'; exit 1 ;;
esac
