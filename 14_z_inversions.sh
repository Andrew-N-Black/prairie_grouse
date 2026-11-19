#!/bin/bash
# =============================================================================
# SLURM JOB SUBMISSION: LARGE INVERSIONS ON THE Z CHROMOSOME
# Step 14 — requires step 12 (pangenome/chr_Z graph) and, for the reference
# path, step 13 (structural_variants/chr_Z/ref.txt).
#
# Step 13's variant calls split bubbles larger than 100 kb, so inversions above
# that size were not resolved. This step finds them by whole-chromosome alignment:
#   1. odgi paths -f: the sequence of every Z path in the graph (PanSN names,
#      the same pieces the graph was built from), used for alignment, gaps
#      (>= 10 N) and path ends
#   2. minimap2 -x asm10: every Z path aligned to the reference Z path (the same
#      LEPC male haplotype as step 13); each alignment block is forward or
#      reverse relative to the reference. (odgi untangle, the graph-native
#      route, aborts on this graph with an sdsl assertion, so alignment is used.)
#   3. optional: the reference Z aligned to the rock ptarmigan Z, so events are
#      also reported in the ptarmigan coordinates the manuscript uses and tested
#      for overlap with the 41-44 Mb segment
#   4. 14_z_inversions.py: blocks running opposite to the rest of their path,
#      merged and >= MIN_INV on the reference, clustered across haplotypes into
#      events and classified as diagnostic / fixed between a species pair /
#      polymorphic / private, with a breakpoint-quality flag
#
# USAGE (from the repo root; put minimap2 on PATH first)
#   export PATH=/home/blackan/.conda/envs/zindel/bin:$PATH
#   PTARMIGAN_FA=/path/to/ptarmigan.fna sbatch --export=ALL 14_z_inversions.sh
#   STAGE=analyze sbatch --export=ALL 14_z_inversions.sh     # rerun step 4 only
#
# OUTPUT: ${OUT}/{summary.txt, inversion_events.tsv, inversion_calls.tsv,
#                 validation_candidates.tsv, z_to_ref.paf, ref_to_ptarmigan.paf}
# Validate any event that separates species with HiFi reads spanning its
# breakpoints and with Hi-C contact maps before reporting it.
# =============================================================================
#SBATCH --job-name=z_inversions
#SBATCH -A dewoody
#SBATCH -p cpu
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=32
#SBATCH --mem=64G
#SBATCH -t 06:00:00
#SBATCH --output=logs/%x_%j.out
#SBATCH --error=logs/%x_%j.err
set -euo pipefail
unset LD_PRELOAD

PROJECT_DIR="${CLUSTER_SCRATCH}/GROUSE/grouse_asm"
GRAPH_DIR="${PROJECT_DIR}/pangenome/chr_Z"
OUT="${OUT:-${PROJECT_DIR}/z_inversions}"
PGGB_SIF="${PGGB_SIF:-${PROJECT_DIR}/containers/pggb.sif}"
REF="${REF:-$(cut -f1 "${PROJECT_DIR}/structural_variants/chr_Z/ref.txt" 2>/dev/null || true)}"
PTARMIGAN_FA="${PTARMIGAN_FA:-}"        # rock ptarmigan genome FASTA (GCF_023343835.1); optional
PTARMIGAN_Z="${PTARMIGAN_Z:-NC_064472.1}"
MIN_BLOCK="${MIN_BLOCK:-5000}"          # ignore alignment blocks shorter than this (bp)
MIN_MAPQ="${MIN_MAPQ:-5}"
MIN_INV="${MIN_INV:-50000}"             # minimum inversion span on the reference (bp)
MERGE="${MERGE:-10000}"                 # merge inverted blocks closer than this (bp)
FLANK="${FLANK:-10000}"                 # breakpoint within this of a gap or path end = flagged
THREADS="${SLURM_CPUS_PER_TASK:-8}"
STAGE="${STAGE:-all}"
REPO="${SLURM_SUBMIT_DIR:-$(pwd)}"

APPTAINER_BIN="$(command -v apptainer || command -v singularity || true)"
cx() { local sif="$1"; shift; "$APPTAINER_BIN" exec --cleanenv --bind "${CLUSTER_SCRATCH}" "$sif" "$@"; }

mkdir -p "$OUT" logs
[[ -n "$REF" ]] || { echo "ERROR: no reference path (set REF or run step 13 first)"; exit 1; }
echo "$REF" > "${OUT}/ref.txt"
echo ">>> Z inversions: reference ${REF}  $(date)"

grab() {   # grab <fasta> <sequence name> <output>: extract one sequence
    python3 - "$1" "$2" "$3" << 'PYEOF'
import sys, gzip
fa, want, dest = sys.argv[1:4]
op = gzip.open if fa.endswith(".gz") else open
keep, n = False, 0
with op(fa, "rt") as fh, open(dest, "w") as o:
    for line in fh:
        if line.startswith(">"):
            keep = line[1:].split()[0] == want
            n += keep
        if keep:
            o.write(line)
if not n:
    sys.exit(f"ERROR: {want} not found in {fa}")
PYEOF
}

if [[ "$STAGE" == "all" ]]; then
    command -v minimap2 >/dev/null || { echo "ERROR: minimap2 not on PATH (export PATH=<conda env>/bin:\$PATH)"; exit 1; }
    OG=$(find "$GRAPH_DIR" -maxdepth 1 -name "*.final.og" | head -n1)
    [[ -s "$OG" ]] || { echo "ERROR: no *.final.og in ${GRAPH_DIR}"; exit 1; }

    # 1. path sequences
    cx "$PGGB_SIF" odgi paths -i "$OG" -f > "${OUT}/z_paths.fa"
    grep -c '>' "${OUT}/z_paths.fa" | xargs echo "  Z paths:"
    grab "${OUT}/z_paths.fa" "$REF" "${OUT}/ref_Z.fa"

    # 2. every path -> reference Z; convert to the table 14_z_inversions.py reads
    minimap2 -x asm10 -t "$THREADS" --secondary=no "${OUT}/ref_Z.fa" "${OUT}/z_paths.fa" > "${OUT}/z_to_ref.paf"
    awk -F'\t' -v r="$REF" -v mb="$MIN_BLOCK" -v mq="$MIN_MAPQ" '
        BEGIN { OFS = "\t"; print "query.name","query.start","query.end","ref.name","ref.start","ref.end","score","inv","self.cov","n.th" }
        $1 != r && $6 == r && $12 >= mq && ($4 - $3) >= mb {
            print $1, $3, $4, $6, $8, $9, $10 / $11, ($5 == "-" ? 1 : 0), 1, 1 }' \
        "${OUT}/z_to_ref.paf" > "${OUT}/untangle.tsv"
    echo "  alignment blocks kept: $(($(wc -l < "${OUT}/untangle.tsv") - 1)) from $(wc -l < "${OUT}/z_to_ref.paf") in the PAF"

    # 3. optional: reference Z -> ptarmigan Z, for manuscript coordinates
    rm -f "${OUT}/ref_to_ptarmigan.paf"
    if [[ -n "$PTARMIGAN_FA" ]]; then
        grab "$PTARMIGAN_FA" "$PTARMIGAN_Z" "${OUT}/ptarmigan_Z.fa"
        minimap2 -x asm10 -t "$THREADS" --secondary=no "${OUT}/ptarmigan_Z.fa" "${OUT}/ref_Z.fa" > "${OUT}/ref_to_ptarmigan.paf"
        rm -f "${OUT}/ptarmigan_Z.fa"
        echo "  ptarmigan liftover: $(wc -l < "${OUT}/ref_to_ptarmigan.paf") alignment blocks"
    else
        echo "  (skipping ptarmigan coordinates: set PTARMIGAN_FA)"
    fi
fi

# 4. call, cluster and classify
python3 "${REPO}/14_z_inversions.py" --dir "$OUT" --min-inv "$MIN_INV" --merge "$MERGE" --flank "$FLANK"
echo ">>> done: $(date)"
