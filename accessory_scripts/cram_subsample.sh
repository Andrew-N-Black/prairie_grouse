#!/bin/bash
#SBATCH --job-name=grouse_subsample
#SBATCH -A dewoody
#SBATCH -t 15:00:00
#SBATCH --nodes=1
#SBATCH --ntasks=4
#SBATCH --mem=40G
#SBATCH -p cpu
#SBATCH --mail-type=BEGIN,END,FAIL
#SBATCH --mail-user=blackan@purdue.edu

ml biocontainers samtools
set -euo pipefail

# ============================================================
# Subsample BAMs to a target mean depth, convert to CRAM,
# then QC every output against the target.
# ============================================================

TARGET_DEPTH=4.66
TOLERANCE=0.15          # flag if realized depth is off target by >15%
REF="/scratch/gautschi/blackan/GROUSE/grouse_asm/ref/GCF_026119805.1_pur_lepc_1.0_genomic.fna"
CRAM_DIR="/scratch/gautschi/blackan/GROUSE/output_shotgun/preprocessing/markduplicates/ALL"
OUT_DIR="/scratch/gautschi/blackan/GROUSE/output_shotgun/preprocessing/markduplicates/subsampled"
SEED=42                 # fixed seed for reproducibility across samples
THREADS=4

mkdir -p "$OUT_DIR"
LOG="$OUT_DIR/subsample_log.tsv"
QC_LOG="$OUT_DIR/qc_report.tsv"
 
echo -e "sample\toriginal_depth\ttarget_depth\tfraction\tseedfrac\tstatus" > "$LOG"
echo -e "sample\trealized_depth\ttarget_depth\tpct_diff\tflag" > "$QC_LOG"
 
mean_depth () {
    # Genome-wide mean depth from samtools coverage (length-weighted across contigs)
    # samtools coverage has no short -T flag; reference is passed via --reference
    samtools coverage --reference "$REF" "$1" | awk 'NR>1 {sum+=$7*($3-$2+1); len+=($3-$2+1)} END {print sum/len}'
}
 
# ---------------- Step 1: subsample + convert to CRAM ----------------
for cram in "$CRAM_DIR"/*md.dedup_q20.cram; do
    sample=$(basename "$cram" .md.dedup_q20.cram)
    echo "=== Processing $sample ==="
 
    current_depth=$(mean_depth "$cram")
    frac=$(awk -v t="$TARGET_DEPTH" -v c="$current_depth" 'BEGIN {f=t/c; if (f>1) f=1; printf "%.4f", f}')
 
    if (( $(awk -v f="$frac" 'BEGIN {print (f>=1)}') )); then
        echo "  WARNING: ${sample} depth (${current_depth}x) already at/below target — copying as-is"
        cp "$cram" "$OUT_DIR/${sample}.subsampled.cram"
        samtools index "$OUT_DIR/${sample}.subsampled.cram"
        echo -e "${sample}\t${current_depth}\t${TARGET_DEPTH}\t1.0000\tN/A\tcopied_no_subsample" >> "$LOG"
    else
        seedfrac=$(awk -v s="$SEED" -v f="$frac" 'BEGIN {printf "%d.%s", s, substr(f,3)}')
        # Read CRAM (-T ref) and write CRAM (-C -T ref) directly — no BAM intermediate needed
        samtools view -@ "$THREADS" -T "$REF" -s "$seedfrac" -C \
            -o "$OUT_DIR/${sample}.subsampled.cram" "$cram"
        samtools index "$OUT_DIR/${sample}.subsampled.cram"
        echo -e "${sample}\t${current_depth}\t${TARGET_DEPTH}\t${frac}\t${seedfrac}\tsubsampled" >> "$LOG"
    fi
done
 
echo ""
echo "Subsampling complete. Log written to $LOG"
 
# ---------------- Step 2: QC pass on outputs ----------------
echo ""
echo "=== Running QC on subsampled outputs ==="
 
for cram in "$OUT_DIR"/*.subsampled.cram; do
    sample=$(basename "$cram" .subsampled.cram)
    realized_depth=$(mean_depth "$cram")
 
    pct_diff=$(awk -v r="$realized_depth" -v t="$TARGET_DEPTH" \
        'BEGIN {printf "%.2f", ((r-t)/t)*100}')
 
    flag="OK"
    if (( $(awk -v p="$pct_diff" -v tol="$TOLERANCE" \
        'BEGIN {print (p<-tol*100 || p>tol*100)}') )); then
        flag="FLAGGED"
    fi
 
    echo -e "${sample}\t${realized_depth}\t${TARGET_DEPTH}\t${pct_diff}%\t${flag}" >> "$QC_LOG"
    echo "  ${sample}: realized=${realized_depth}x (${pct_diff}% from target) [${flag}]"
done
 
echo ""
echo "QC complete. Report written to $QC_LOG"
n_flagged=$(awk -F'\t' 'NR>1 && $5=="FLAGGED"' "$QC_LOG" | wc -l)
echo "Samples flagged as off-target (>${TOLERANCE}x tolerance): $n_flagged"
if [ "$n_flagged" -gt 0 ]; then
    echo "  Review $QC_LOG and consider re-running flagged samples with an adjusted fraction."
fi
