#!/bin/bash
# =============================================================================
# SLURM ARRAY JOB: COMBINED-COHORT ALIGNMENT TO LEPC REFERENCE
# Step 05 of the replication/extension of Black et al. (PNAS Nexus)
# LEPC-popgen analysis (https://github.com/Andrew-N-Black/LEPC-popgen).
#
# Combines two cohorts into one alignment array:
#   - OLD: n=468 raw paired-end FASTQs from PRJNA986511 (the original PNAS
#     Nexus samples), auto-discovered from PRJNA986511_RAW_DIR.
#   - NEW: additional short-read population samples from
#     shortread_manifest.tsv (built for 04_shortread_introgression.sh).
#
# Per sample: bwa mem -> sort -> GATK4 MarkDuplicates -> quality filter ->
# per-sample mean depth. This deliberately DIFFERS from the original
# alignment.sh in two ways (per explicit instruction, not an oversight):
#   - GATK4 MarkDuplicates instead of GATK 3.6.0 RealignerTargetCreator/
#     IndelRealigner — GATK 3.6.0 is long deprecated and unavailable here;
#     local realignment is considered obsolete with modern callers anyway.
#   - No mappability-mask / RepeatMasker BED restriction — filtering is
#     quality-based only (matches the original's samtools view flags).
#
# Depth here is NOT yet the final analysis depth — 06_downsample_and_
# finalize.sh compares depths across both cohorts and downsamples any NEW
# sample that exceeds the OLD cohort's empirical mean before anything
# downstream (heterozygosity/PCA/ROH) runs, so both cohorts are analyzed
# at matched coverage.
#
# USAGE:
#   1. Build the combined sample list + prep the reference ONCE:
#        sbatch --array=0-0 05_combined_alignment_array.sh --build-list-only
#
#   2. Submit the full array (the prep step above prints the exact N):
#        N=$(grep -vc '^#' <PROJECT_DIR>/combined_samples.tsv)
#        sbatch --array=0-$((N-1))%20 05_combined_alignment_array.sh
# =============================================================================
#SBATCH --job-name=lepc_align
#SBATCH --output=logs/%x_%A_%a.out
#SBATCH --error=logs/%x_%A_%a.err
#SBATCH -A fnrdewoody
#SBATCH -t 2-00:00:00
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --mem=32G
#SBATCH -p cpu
#SBATCH --mail-type=BEGIN,END,FAIL
#SBATCH --mail-user=${USER}@purdue.edu

# =============================================================================
# ENVIRONMENT SETUP
# =============================================================================
set -euo pipefail

ml biocontainers
ml bwa
ml samtools/1.22.1
ml gatk4

# =============================================================================
# USER-DEFINED VARIABLES
# =============================================================================
PROJECT_DIR="${CLUSTER_SCRATCH}/LEPC"
REF_DIR="${PROJECT_DIR}/ref"

# LEPC (Tympanuchus pallidicinctus) RefSeq reference — GCF_026119805.1 /
# pur_lepc_1.0. Downloaded by 01_download_reference_genomes.sh; this script
# only builds the tool-specific indices it needs (bwa/gatk), not the raw
# download itself.
REF_FASTA="${REF_DIR}/GCF_026119805.1_pur_lepc_1.0_genomic.fa"

PRJNA986511_RAW_DIR="${CLUSTER_SCRATCH}/PRJNA986511/raw"
SHORTREAD_MANIFEST="${SLURM_SUBMIT_DIR}/shortread_manifest.tsv"

COMBINED_SAMPLES="${PROJECT_DIR}/combined_samples.tsv"

ALIGN_DIR="${PROJECT_DIR}/aligned"
DEPTH_DIR="${PROJECT_DIR}/depth"

THREADS=$SLURM_CPUS_PER_TASK

mkdir -p logs "$REF_DIR" "$ALIGN_DIR" "$DEPTH_DIR"

# =============================================================================
# BLOCK A: ONE-TIME PREP — reference prep + combined sample list.
# Run once before the array (see USAGE above). Skipped automatically on
# subsequent runs if outputs already exist.
# =============================================================================
prep_reference() {
    if [[ ! -f "$REF_FASTA" ]]; then
        echo "ERROR: LEPC reference not found: ${REF_FASTA}"
        echo "Run 01_download_reference_genomes.sh first."
        exit 1
    fi
    if [[ ! -f "${REF_FASTA}.fai" ]]; then
        echo ">>> Indexing reference (samtools faidx)"
        samtools faidx "$REF_FASTA"
    fi
    if [[ ! -f "${REF_FASTA}.bwt" ]]; then
        echo ">>> Indexing reference (bwa index)"
        bwa index "$REF_FASTA"
    fi
    local dict="${REF_FASTA%.fa}.dict"
    if [[ ! -f "$dict" ]]; then
        echo ">>> Building sequence dictionary (gatk CreateSequenceDictionary)"
        gatk CreateSequenceDictionary -R "$REF_FASTA" -O "$dict"
    fi
}

build_combined_samples() {
    echo ">>> Building combined sample list -> ${COMBINED_SAMPLES}"
    echo -e "sample_id\tcohort\tfastq_r1\tfastq_r2" > "$COMBINED_SAMPLES"

    # --- OLD cohort: auto-discover paired FASTQs from PRJNA986511 ---
    if [[ ! -d "$PRJNA986511_RAW_DIR" ]]; then
        echo "ERROR: PRJNA986511 raw directory not found: ${PRJNA986511_RAW_DIR}"
        exit 1
    fi

    local n_old=0
    while IFS= read -r r1; do
        local basename sample r2
        basename=$(basename "$r1")
        case "$basename" in
            *_R1_001.fastq.gz) sample="${basename%_R1_001.fastq.gz}"; r2="${r1/_R1_001.fastq.gz/_R2_001.fastq.gz}" ;;
            *_R1.fastq.gz)     sample="${basename%_R1.fastq.gz}";     r2="${r1/_R1.fastq.gz/_R2.fastq.gz}" ;;
            *_1.fastq.gz)      sample="${basename%_1.fastq.gz}";      r2="${r1/_1.fastq.gz/_2.fastq.gz}" ;;
            *) echo "  WARNING: unrecognized naming pattern for ${basename} — skipping" >&2; continue ;;
        esac
        if [[ ! -f "$r2" ]]; then
            echo "  WARNING: mate not found for ${r1} — skipping" >&2
            continue
        fi
        sample=$(echo "$sample" | tr -c 'A-Za-z0-9_-' '_')
        echo -e "${sample}\tOLD\t${r1}\t${r2}" >> "$COMBINED_SAMPLES"
        n_old=$((n_old + 1))
    done < <(find "$PRJNA986511_RAW_DIR" -maxdepth 1 -type f \( -name "*_R1_001.fastq.gz" -o -name "*_R1.fastq.gz" -o -name "*_1.fastq.gz" \) | sort)

    echo "  OLD cohort: ${n_old} sample pairs"

    # --- NEW cohort: from shortread_manifest.tsv ---
    if [[ ! -f "$SHORTREAD_MANIFEST" ]]; then
        echo "ERROR: Short-read manifest not found: ${SHORTREAD_MANIFEST}"
        exit 1
    fi
    local n_new=0
    while IFS=$'\t' read -r sample species r1 r2; do
        [[ -z "$sample" ]] && continue
        echo -e "${sample}\tNEW\t${r1}\t${r2}" >> "$COMBINED_SAMPLES"
        n_new=$((n_new + 1))
    done < <(grep -v '^#' "$SHORTREAD_MANIFEST" | tail -n +2)

    echo "  NEW cohort: ${n_new} sample pairs"

    local total=$((n_old + n_new))
    echo "  Total: ${total} samples -> ${COMBINED_SAMPLES}"
    echo ""
    echo ">>> Submit the full array with:"
    echo "    sbatch --array=0-$((total - 1))%20 05_combined_alignment_array.sh"
}

if [[ "${1:-}" == "--build-list-only" ]]; then
    prep_reference
    build_combined_samples
    echo ">>> Prep finished."
    exit 0
fi

# Safety net (same caveat as elsewhere in this pipeline: run --build-list-only
# first to avoid every array task racing on this simultaneously).
prep_reference

if [[ ! -f "$COMBINED_SAMPLES" ]]; then
    echo "ERROR: ${COMBINED_SAMPLES} not found."
    echo "Run: sbatch --array=0-0 05_combined_alignment_array.sh --build-list-only"
    exit 1
fi

# =============================================================================
# RESOLVE SAMPLE FOR THIS ARRAY TASK
# =============================================================================
if [[ -z "${SLURM_ARRAY_TASK_ID:-}" ]]; then
    echo "ERROR: SLURM_ARRAY_TASK_ID is not set."
    echo "Submit with: sbatch --array=0-N 05_combined_alignment_array.sh"
    exit 1
fi

mapfile -t ROWS < <(grep -v '^#' "$COMBINED_SAMPLES" | tail -n +2)
LINE="${ROWS[$SLURM_ARRAY_TASK_ID]:-}"

if [[ -z "$LINE" ]]; then
    echo "ERROR: No sample at index ${SLURM_ARRAY_TASK_ID} in ${COMBINED_SAMPLES}"
    exit 1
fi

IFS=$'\t' read -r SAMPLE COHORT R1 R2 <<< "$LINE"

echo ">>> Array task ${SLURM_ARRAY_TASK_ID} -> sample: ${SAMPLE} (${COHORT})"
echo ">>> Reads: ${R1} / ${R2}"
echo ">>> CPUs : ${THREADS}"
echo ">>> Start: $(date)"

for f in "$R1" "$R2"; do
    if [[ ! -f "$f" ]]; then
        echo "ERROR: Input file not found: ${f}"
        exit 1
    fi
done

# =============================================================================
# STEP 1: bwa mem alignment + sort
# =============================================================================
echo ">>> Step 1: bwa mem alignment"

SORTED_BAM="${ALIGN_DIR}/${SAMPLE}.sorted.bam"
RG="@RG\tID:${SAMPLE}\tSM:${SAMPLE}\tPL:ILLUMINA\tLB:${SAMPLE}"

bwa mem -t "$THREADS" -M -R "$RG" "$REF_FASTA" "$R1" "$R2" \
    | samtools sort -@ "$THREADS" -o "$SORTED_BAM" -
samtools index "$SORTED_BAM"

# =============================================================================
# STEP 2: GATK4 MarkDuplicates
# =============================================================================
echo ">>> Step 2: GATK4 MarkDuplicates"

DEDUP_BAM="${ALIGN_DIR}/${SAMPLE}.dedup.bam"
gatk MarkDuplicates \
    -I "$SORTED_BAM" \
    -O "$DEDUP_BAM" \
    -M "${ALIGN_DIR}/${SAMPLE}.dedup_metrics.txt" \
    --OPTICAL_DUPLICATE_PIXEL_DISTANCE 2500
samtools index "$DEDUP_BAM"

# =============================================================================
# STEP 3: Quality filter (matches the original alignment.sh's samtools view
# flags exactly — MAPQ>=30, proper pairs only, exclude unmapped/secondary/
# qcfail/dup/supplementary; no BED restriction, per explicit instruction)
# =============================================================================
echo ">>> Step 3: Quality filter (MAPQ>=30, -F 3844 -f 2)"

FILT_BAM="${ALIGN_DIR}/${SAMPLE}_filt.bam"
samtools view -@ "$THREADS" -q 30 -b -F 3844 -f 2 "$DEDUP_BAM" > "$FILT_BAM"
samtools index "$FILT_BAM"

# =============================================================================
# STEP 4: Per-sample depth (mean depth + breadth >=1x), matching the
# original script's coverage-stats calculation
# =============================================================================
echo ">>> Step 4: Coverage statistics"

DEPTH_FILE="${DEPTH_DIR}/${SAMPLE}.depth.txt"
samtools depth -a "$FILT_BAM" | awk -v sample="$SAMPLE" -v cohort="$COHORT" '
    { sum += $3; n++; if ($3 >= 1) covered++ }
    END {
        mean = (n > 0) ? sum / n : 0
        breadth = (n > 0) ? 100 * covered / n : 0
        printf "%s\t%s\t%.4f\t%.2f\n", sample, cohort, mean, breadth
    }
' > "$DEPTH_FILE"

echo "  $(cat "$DEPTH_FILE")"

# Clean up intermediate (pre-dedup, pre-filter) BAMs — the filtered BAM is
# what everything downstream uses; keeping the sorted+dedup intermediates
# for ~500 samples would meaningfully add to scratch usage for no benefit.
rm -f "$SORTED_BAM" "${SORTED_BAM}.bai"

echo ">>> Sample ${SAMPLE} (${COHORT}) complete — $(date)"
