#!/bin/bash
# =============================================================================
# SLURM JOB SUBMISSION: DEPTH-MATCH DOWNSAMPLING + FINAL BAM LIST
# Step 06 — requires 05_combined_alignment_array.sh to have completed for
# every sample (reads its per-sample depth files).
#
# The OLD cohort (n=468, PRJNA986511) is low coverage (originally reported
# ~6x). The NEW cohort's actual depth is whatever it turns out to be once
# sequenced — likely higher. Mixing cohorts at different depths would bias
# heterozygosity/PCA/ROH (all coverage-sensitive), so this script:
#   1. Computes the OLD cohort's TARGET_DEPTH empirically from its own
#      measured depths (not a hardcoded "6" — more rigorous and
#      self-consistent with whatever the actual data show).
#   2. Downsamples only NEW-cohort BAMs that exceed TARGET_DEPTH, via
#      `samtools view -s <seed>.<fraction>`.
#   3. Leaves OLD-cohort BAMs and any NEW sample already at/below
#      TARGET_DEPTH untouched.
#   4. Writes final_bamlist.txt — the single input every downstream script
#      (07, 08) consumes.
#
# USAGE:
#   sbatch 06_downsample_and_finalize.sh
# =============================================================================
#SBATCH --job-name=lepc_downsample
#SBATCH --output=logs/%x_%j.out
#SBATCH --error=logs/%x_%j.err
#SBATCH -A fnrdewoody
#SBATCH -t 1-00:00:00
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=16
#SBATCH --mem=32G
#SBATCH -p cpu
#SBATCH --mail-type=BEGIN,END,FAIL
#SBATCH --mail-user=${USER}@purdue.edu

# =============================================================================
# ENVIRONMENT SETUP
# =============================================================================
set -euo pipefail

ml biocontainers
ml samtools/1.22.1

# =============================================================================
# USER-DEFINED VARIABLES
# =============================================================================
PROJECT_DIR="${CLUSTER_SCRATCH}/LEPC"
ALIGN_DIR="${PROJECT_DIR}/aligned"
DEPTH_DIR="${PROJECT_DIR}/depth"
DOWNSAMPLE_DIR="${PROJECT_DIR}/downsampled"

FINAL_BAMLIST="${PROJECT_DIR}/final_bamlist.txt"
DEPTH_SUMMARY="${PROJECT_DIR}/depth_summary.tsv"

# Fixed seed for reproducible downsampling (samtools view -s SEED.FRACTION)
DOWNSAMPLE_SEED=42

THREADS=$SLURM_CPUS_PER_TASK

mkdir -p logs "$DOWNSAMPLE_DIR"

echo ">>> 06_downsample_and_finalize.sh"
echo ">>> Start time: $(date)"

# =============================================================================
# STEP 1: Aggregate per-sample depth files
# =============================================================================
echo ">>> Step 1: Aggregating depth files from ${DEPTH_DIR}"

N_DEPTH_FILES=$(find "$DEPTH_DIR" -name "*.depth.txt" | wc -l)
if [[ "$N_DEPTH_FILES" -eq 0 ]]; then
    echo "ERROR: No depth files found in ${DEPTH_DIR}"
    echo "Run 05_combined_alignment_array.sh for all samples first."
    exit 1
fi

echo -e "sample_id\tcohort\tmean_depth\tbreadth_pct" > "$DEPTH_SUMMARY"
find "$DEPTH_DIR" -name "*.depth.txt" -exec cat {} + >> "$DEPTH_SUMMARY"

N_SAMPLES=$((N_DEPTH_FILES))
echo "  ${N_SAMPLES} samples with depth data"

# =============================================================================
# STEP 2: Compute TARGET_DEPTH from the OLD cohort's empirical mean depth
# =============================================================================
TARGET_DEPTH=$(awk -F'\t' 'NR>1 && $2=="OLD" { sum+=$3; n++ } END { if (n>0) printf "%.4f", sum/n; else print "0" }' "$DEPTH_SUMMARY")
N_OLD=$(awk -F'\t' 'NR>1 && $2=="OLD"' "$DEPTH_SUMMARY" | wc -l)
N_NEW=$(awk -F'\t' 'NR>1 && $2=="NEW"' "$DEPTH_SUMMARY" | wc -l)

echo ">>> Step 2: TARGET_DEPTH = ${TARGET_DEPTH}x (empirical mean of ${N_OLD} OLD-cohort samples)"
echo "    Cross-check: the original paper reports ~6x mean depth for this cohort —"
echo "    if TARGET_DEPTH above is wildly different from 6, investigate before proceeding."

if [[ "$N_OLD" -eq 0 ]]; then
    echo "ERROR: No OLD-cohort samples found — cannot compute a target depth."
    exit 1
fi
if [[ "$N_NEW" -eq 0 ]]; then
    echo "  NOTE: no NEW-cohort samples found yet (expected while sequencing is pending)."
    echo "  This script will still finalize a bamlist from the OLD cohort alone."
fi

# =============================================================================
# STEP 3: Downsample NEW-cohort BAMs exceeding TARGET_DEPTH; pass everything
# else through unchanged.
# =============================================================================
echo ">>> Step 3: Downsampling NEW-cohort samples exceeding target depth"

DOWNSAMPLE_LOG="${PROJECT_DIR}/downsample_log.tsv"
echo -e "sample_id\tcohort\tmeasured_depth\ttarget_depth\tfraction_kept\taction" > "$DOWNSAMPLE_LOG"

> "$FINAL_BAMLIST"

while IFS=$'\t' read -r sample cohort depth breadth; do
    [[ "$sample" == "sample_id" ]] && continue
    filt_bam="${ALIGN_DIR}/${sample}_filt.bam"

    if [[ ! -f "$filt_bam" ]]; then
        echo "  WARNING: expected BAM not found for ${sample} — skipping: ${filt_bam}" >&2
        continue
    fi

    if [[ "$cohort" == "OLD" ]]; then
        echo -e "${sample}\tOLD\t${depth}\t${TARGET_DEPTH}\t1.0000\tunchanged" >> "$DOWNSAMPLE_LOG"
        echo "$filt_bam" >> "$FINAL_BAMLIST"
        continue
    fi

    # NEW cohort: downsample only if measured depth exceeds target
    exceeds=$(awk -v d="$depth" -v t="$TARGET_DEPTH" 'BEGIN { print (d > t) ? 1 : 0 }')
    if [[ "$exceeds" -eq 1 ]]; then
        fraction=$(awk -v d="$depth" -v t="$TARGET_DEPTH" 'BEGIN { printf "%.4f", t/d }')
        # Build the full "SEED.FRACTION" string in one awk call rather than
        # stripping a leading "0." in bash — that string trick silently
        # produces a malformed samtools argument (e.g. "42.1.0000") when the
        # fraction rounds up to 1.0000 (depth only marginally above target).
        seedfrac=$(awk -v d="$depth" -v t="$TARGET_DEPTH" -v seed="$DOWNSAMPLE_SEED" \
            'BEGIN { f = t/d; if (f >= 1) f = 0.9999; printf "%d.%04d", seed, int(f*10000) }')
        ds_bam="${DOWNSAMPLE_DIR}/${sample}_ds.bam"
        if [[ ! -f "$ds_bam" ]]; then
            samtools view -@ "$THREADS" -s "$seedfrac" -b "$filt_bam" > "$ds_bam"
            samtools index "$ds_bam"
        fi
        echo -e "${sample}\tNEW\t${depth}\t${TARGET_DEPTH}\t${fraction}\tdownsampled" >> "$DOWNSAMPLE_LOG"
        echo "$ds_bam" >> "$FINAL_BAMLIST"
    else
        echo -e "${sample}\tNEW\t${depth}\t${TARGET_DEPTH}\t1.0000\tunchanged (already <= target)" >> "$DOWNSAMPLE_LOG"
        echo "$filt_bam" >> "$FINAL_BAMLIST"
    fi
done < "$DEPTH_SUMMARY"

N_FINAL=$(wc -l < "$FINAL_BAMLIST")

echo ""
echo ">>> Downsampling complete."
echo "    Target depth       : ${TARGET_DEPTH}x"
echo "    Samples in final list: ${N_FINAL}"
echo "    Downsample log      : ${DOWNSAMPLE_LOG}"
echo "    Final BAM list      : ${FINAL_BAMLIST}"
echo ""
echo ">>> Next: submit 07_heterozygosity_array.sh with:"
echo "    sbatch --array=0-$((N_FINAL - 1))%20 07_heterozygosity_array.sh"
echo ">>> End time: $(date)"
