#!/bin/bash
# =============================================================================
# SLURM JOB SUBMISSION: LARGE INVERSIONS ON THE Z CHROMOSOME
# Step 06 — requires step 03 (pangenome/chr_Z graph) and, for the reference
# path, step 05 (structural_variants/chr_Z/ref.txt).
#
# Step 05's variant calls split bubbles larger than 100 kb, so inversions above
# that size were not resolved. This step reads them directly from the Z graph:
#   1. odgi untangle: every Z path is segmented against one reference Z path
#      (the same LEPC male haplotype as step 05), each segment forward or
#      reverse relative to it
#   2. odgi paths -f: every Z path's sequence, to locate gaps (>= 10 N) and
#      path ends, where a mis-oriented scaffold would fake an inversion
#   3. optional: the reference Z aligned to the rock ptarmigan Z (minimap2), so
#      events are also reported in the ptarmigan coordinates the manuscript
#      uses, and tested for overlap with the 41-44 Mb segment
#   4. 06_z_inversions.py: inverted blocks >= MIN_INV per haplotype, clustered
#      across haplotypes into events, classified as diagnostic / fixed between a
#      species pair / polymorphic / private, with a breakpoint-quality flag
#
# USAGE (from the repo root)
#   sbatch 06_z_inversions.sh                         # everything
#   STAGE=analyze sbatch --export=ALL 06_z_inversions.sh   # rerun step 4 only (new thresholds)
#   Faster on the full Z: one array task per query path, then the rest:
#     bash 06_z_inversions.sh split        # submits the array + a dependent STAGE=post job
#   (UNTANGLE_M sets odgi untangle -m; default 256, the split mode uses 1000)
#   Put minimap2 on PATH first for the ptarmigan coordinates, e.g.
#     export PATH=/home/blackan/.conda/envs/zindel/bin:$PATH
#
# OUTPUT: ${OUT}/{summary.txt, inversion_events.tsv, inversion_calls.tsv,
#                 validation_candidates.tsv, untangle.tsv, ref_to_ptarmigan.paf}
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
#SBATCH -t 12:00:00
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
MIN_INV="${MIN_INV:-50000}"             # minimum inversion span on the reference (bp)
MERGE="${MERGE:-10000}"                 # merge inverted segments closer than this (bp)
FLANK="${FLANK:-10000}"                 # breakpoint within this of a gap or path end = flagged
THREADS="${SLURM_CPUS_PER_TASK:-8}"
STAGE="${STAGE:-${1:-all}}"
UNTANGLE_M="${UNTANGLE_M:-256}"
REPO="${SLURM_SUBMIT_DIR:-$(pwd)}"

APPTAINER_BIN="$(command -v apptainer || command -v singularity || true)"
cx() { local sif="$1"; shift; "$APPTAINER_BIN" exec --cleanenv --bind "${CLUSTER_SCRATCH}" "$sif" "$@"; }

mkdir -p "$OUT" logs
[[ -n "$REF" ]] || { echo "ERROR: no reference path (set REF or run step 05 first)"; exit 1; }
echo "$REF" > "${OUT}/ref.txt"
echo ">>> Z inversions: reference ${REF}  $(date)"

OG=$(find "$GRAPH_DIR" -maxdepth 1 -name "*.final.og" | head -n1)
[[ -s "$OG" ]] || { echo "ERROR: no *.final.og in ${GRAPH_DIR}"; exit 1; }
list_paths() {
    cx "$PGGB_SIF" odgi paths -i "$OG" -L > "${OUT}/paths.txt"
    grep -q -x -F "$REF" "${OUT}/paths.txt" || { echo "ERROR: ${REF} is not a path in ${OG}"; exit 1; }
    echo "$REF" > "${OUT}/target.txt"
    grep -v -x -F "$REF" "${OUT}/paths.txt" > "${OUT}/queries.txt"
}

# split mode (login node): one array task per query path, then STAGE=post
if [[ "$STAGE" == "split" ]]; then
    list_paths
    N=$(wc -l < "${OUT}/queries.txt"); mkdir -p "${OUT}/untangle_parts"
    jid=$(sbatch --parsable --export=ALL,STAGE=one,OUT="$OUT",REF="$REF",UNTANGLE_M="${UNTANGLE_M:-1000}" \
          --array=0-$((N - 1)) -c 8 --mem=32G -t 12:00:00 --job-name=z_untangle "$0")
    pid=$(sbatch --parsable --export=ALL,STAGE=post,OUT="$OUT",REF="$REF" --dependency=afterok:"$jid" "$0")
    echo "Submitted ${N} untangle tasks (array ${jid}) and the post step (${pid})"
    exit 0
fi

# one array task: untangle a single query path
if [[ "$STAGE" == "one" ]]; then
    Q=$(sed -n "$((SLURM_ARRAY_TASK_ID + 1))p" "${OUT}/queries.txt")
    echo "  untangle ${Q} (-m ${UNTANGLE_M})"
    cx "$PGGB_SIF" odgi untangle -i "$OG" -r "$REF" -q "$Q" -m "$UNTANGLE_M" -j 0 -n 1 -t "$THREADS" \
        > "${OUT}/untangle_parts/part_$(printf %03d "$SLURM_ARRAY_TASK_ID").tsv"
    exit 0
fi

if [[ "$STAGE" == "all" ]]; then
    # 1. untangle every path against the reference path
    list_paths
    echo "  odgi untangle: $(wc -l < "${OUT}/queries.txt") query paths"
    cx "$PGGB_SIF" odgi untangle -i "$OG" -R "${OUT}/target.txt" -Q "${OUT}/queries.txt" \
        -m "$UNTANGLE_M" -j 0 -n 1 -t "$THREADS" > "${OUT}/untangle.tsv"
    echo "  untangle rows: $(grep -vc 'query.name' "${OUT}/untangle.tsv")"
fi

if [[ "$STAGE" == "post" ]]; then
    # join the per-path parts: one header, then every part's rows
    { grep -h -m1 'query.name' "${OUT}"/untangle_parts/part_*.tsv | head -1
      for f in "${OUT}"/untangle_parts/part_*.tsv; do grep -v 'query.name' "$f" || true; done
    } > "${OUT}/untangle.tsv"
    echo "  untangle rows: $(grep -vc 'query.name' "${OUT}/untangle.tsv") from $(ls "${OUT}"/untangle_parts/part_*.tsv | wc -l) parts"
fi

if [[ "$STAGE" == "all" || "$STAGE" == "post" ]]; then

    # 2. path sequences (gaps and lengths)
    cx "$PGGB_SIF" odgi paths -i "$OG" -f > "${OUT}/z_paths.fa"

    # 3. optional: reference Z -> ptarmigan Z, for manuscript coordinates
    rm -f "${OUT}/ref_to_ptarmigan.paf"
    if [[ -n "$PTARMIGAN_FA" ]] && command -v minimap2 >/dev/null; then
        python3 - "$PTARMIGAN_FA" "$PTARMIGAN_Z" "${OUT}/z_paths.fa" "$REF" "$OUT" << 'PYEOF'
import sys, gzip
fa, zname, paths, ref, out = sys.argv[1:6]
def grab(path, want, dest):
    op = gzip.open if path.endswith(".gz") else open
    keep, n = False, 0
    with op(path, "rt") as fh, open(dest, "w") as o:
        for line in fh:
            if line.startswith(">"):
                keep = line[1:].split()[0] == want
                n += keep
            if keep:
                o.write(line)
    if not n:
        sys.exit(f"ERROR: {want} not found in {path}")
grab(fa, zname, f"{out}/ptarmigan_Z.fa")
grab(paths, ref, f"{out}/ref_Z.fa")
PYEOF
        minimap2 -x asm5 -t "$THREADS" "${OUT}/ptarmigan_Z.fa" "${OUT}/ref_Z.fa" > "${OUT}/ref_to_ptarmigan.paf"
        rm -f "${OUT}/ptarmigan_Z.fa" "${OUT}/ref_Z.fa"
        echo "  ptarmigan liftover: $(wc -l < "${OUT}/ref_to_ptarmigan.paf") alignment blocks"
    else
        echo "  (skipping ptarmigan coordinates: set PTARMIGAN_FA and put minimap2 on PATH)"
    fi
fi

# 4. call, cluster and classify
python3 "${REPO}/06_z_inversions.py" --dir "$OUT" --min-inv "$MIN_INV" --merge "$MERGE" --flank "$FLANK"
echo ">>> done: $(date)"
