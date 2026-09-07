#!/bin/bash
# =============================================================================
# SLURM JOB SUBMISSION: nf-core/sarek — LEPC "NEW" vs "OLD" COHORT COMPARISON
# Compares 10 "New" and 10 "Old" LEPC whole-genome paired-end samples,
# both ~20x depth of coverage, via joint germline variant calling.
#
# 20x is comfortably in GATK HaplotypeCaller's reliable range for individual
# hard genotype calls — unlike the ~6x PNAS Nexus cohort in the 05-08
# scripts, which needed ANGSD's genotype-likelihood framework instead. No
# depth-matching/downsampling step is needed here since both batches are
# already comparable depth.
#
# This job only orchestrates Nextflow; the heavy per-sample work is farmed
# out as individual SLURM jobs via the 'slurm' executor profile already
# defined in nextflow.config (reused as-is — it's the same config written
# for nf-core/sarek v3.9.0 on this cluster/account), so the head job itself
# needs only modest resources.
#
# Reuses the LEPC reference already downloaded by
# 01_download_reference_genomes.sh — does not re-fetch it.
#
# BATCH IDENTITY: sarek's samplesheet has no dedicated "batch" column, so
# batch identity is encoded directly into the sample ID (NEW_<name> /
# OLD_<name>). This makes it visible in the joint VCF's sample columns and
# is what Step 3 below uses to split per-batch comparison stats.
#
# INPUT: point NEW_RAW_DIR / OLD_RAW_DIR (below) at your actual raw FASTQ
# directories before submitting — defaults are a guess at this project's
# directory convention, not verified paths.
#
# USAGE:
#   sbatch 09_sarek_new_vs_old.sh
# =============================================================================
#SBATCH --job-name=lepc_new_vs_old
#SBATCH --output=logs/%x_%j.out
#SBATCH --error=logs/%x_%j.err
#SBATCH -A fnrdewoody
#SBATCH -t 14-00:00:00
#SBATCH -p cpu
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=2
#SBATCH --mem=8G
#SBATCH --mail-type=BEGIN,END,FAIL
#SBATCH --mail-user=${USER}@purdue.edu

# =============================================================================
# ENVIRONMENT SETUP
# =============================================================================
set -euo pipefail
mkdir -p logs

# =============================================================================
# USER-DEFINED VARIABLES
# =============================================================================
PROJECT_DIR="${CLUSTER_SCRATCH}/LEPC"

# Reference genome — downloaded by 01_download_reference_genomes.sh, not
# re-fetched here. Same file 05-08 already use.
REF_FASTA="${PROJECT_DIR}/ref/GCF_026119805.1_pur_lepc_1.0_genomic.fa"

# Raw paired-end FASTQ directories — ONE FLAT DIRECTORY per batch, 10
# samples each. Verify/adjust these paths before submitting.
NEW_RAW_DIR="${PROJECT_DIR}/raw_new"
OLD_RAW_DIR="${PROJECT_DIR}/raw_old"

SAMPLESHEET="${PROJECT_DIR}/samplesheet_new_vs_old.csv"
OUTDIR="${PROJECT_DIR}/nf-out"
COMPARE_DIR="${PROJECT_DIR}/new_vs_old"

SAREK_VERSION="3.9.0"

echo ">>> 09_sarek_new_vs_old.sh — LEPC New vs Old cohort comparison"
echo ">>> Reference : ${REF_FASTA}"
echo ">>> New reads : ${NEW_RAW_DIR}"
echo ">>> Old reads : ${OLD_RAW_DIR}"
echo ">>> Output dir: ${OUTDIR}"
echo ">>> Node      : $(hostname)"
echo ">>> Start time: $(date)"

mkdir -p "${OUTDIR}" "${COMPARE_DIR}"

if [[ ! -f "$REF_FASTA" ]]; then
    echo "ERROR: LEPC reference not found: ${REF_FASTA}"
    echo "Run 01_download_reference_genomes.sh first."
    exit 1
fi

# =============================================================================
# STEP 1: Build the sarek samplesheet from both raw directories.
# Handles *_R1_001/_R2_001, *_R1/_R2, and *_1/_2 naming conventions (same
# detection logic used elsewhere in this project). Each pair becomes one
# germline sample (patient == sample, lane 1), prefixed with its batch.
# =============================================================================
echo ">>> Step 1: Building samplesheet from ${NEW_RAW_DIR} and ${OLD_RAW_DIR}"

echo "patient,sample,lane,fastq_1,fastq_2" > "${SAMPLESHEET}"

build_batch() {
    local raw_dir="$1"
    local batch_prefix="$2"
    local n=0

    if [[ ! -d "$raw_dir" ]]; then
        echo "ERROR: ${batch_prefix} raw directory not found: ${raw_dir}"
        exit 1
    fi

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
        local id="${batch_prefix}_${sample}"
        echo "${id},${id},1,${r1},${r2}" >> "${SAMPLESHEET}"
        n=$((n + 1))
    done < <(find "$raw_dir" -maxdepth 1 -type f \( -name "*_R1_001.fastq.gz" -o -name "*_R1.fastq.gz" -o -name "*_1.fastq.gz" \) | sort)

    echo "  ${batch_prefix}: ${n} sample pairs"
    if [[ "$n" -ne 10 ]]; then
        echo "  WARNING: expected 10 ${batch_prefix} samples, found ${n} — check ${raw_dir}" >&2
    fi
}

build_batch "$NEW_RAW_DIR" "NEW"
build_batch "$OLD_RAW_DIR" "OLD"

N_TOTAL=$(($(wc -l < "${SAMPLESHEET}") - 1))
echo "  Wrote ${N_TOTAL} total samples to ${SAMPLESHEET}"

if [[ "${N_TOTAL}" -eq 0 ]]; then
    echo "ERROR: no paired-end FASTQs found in either batch directory" >&2
    exit 1
fi

# =============================================================================
# STEP 2: nf-core/sarek
# fastp adapter/quality clipping -> BWA-MEM2 alignment -> GATK4
# MarkDuplicates + BQSR -> GATK4 HaplotypeCaller (per-sample GVCF) -> joint
# genotyping across all 20 samples. No annotation step — this run is scoped
# to comparing the two batches via genotype calls, not functional annotation.
# =============================================================================
echo ">>> Step 2: nf-core/sarek v${SAREK_VERSION}"

ml biocontainers
ml nf-core/2.11.1

cd "${PROJECT_DIR}"

nextflow run nf-core/sarek \
    -r "${SAREK_VERSION}" \
    -c "${PROJECT_DIR}/nextflow.config" \
    -profile singularity,slurm \
    --input "${SAMPLESHEET}" \
    --outdir "${OUTDIR}" \
    --fasta "${REF_FASTA}" \
    --trim_fastq \
    --tools haplotypecaller \
    --joint_germline \
    --save_reference \
    --email "${USER}@purdue.edu"

# =============================================================================
# STEP 3: Stage the final joint VCF, then split per-batch comparison stats
# (bcftools stats) using the NEW_/OLD_ prefixes embedded in Step 1 — a
# direct, lightweight answer to "how do these two batches compare" without
# a bespoke statistical framework. For anything beyond summary counts
# (PCA, relatedness, etc.), ask for that as a specific follow-up once
# you've seen these numbers.
# =============================================================================
echo ">>> Step 3: Locating final joint VCF + per-batch comparison"

ml bcftools

JOINT_VCF=$(find "${OUTDIR}/variant_calling" -ipath "*joint_variant_calling*" -name "*.vcf.gz" 2>/dev/null | head -n1)

if [[ -z "${JOINT_VCF}" ]]; then
    echo "WARNING: could not auto-locate the joint VCF — check ${OUTDIR}/variant_calling manually"
    echo ">>> End time: $(date)"
    exit 0
fi

cp "${JOINT_VCF}" "${COMPARE_DIR}/LEPC_new_vs_old_joint.vcf.gz"
tabix -f -p vcf "${COMPARE_DIR}/LEPC_new_vs_old_joint.vcf.gz" 2>/dev/null || true
echo "  Joint VCF: ${COMPARE_DIR}/LEPC_new_vs_old_joint.vcf.gz"

bcftools query -l "${COMPARE_DIR}/LEPC_new_vs_old_joint.vcf.gz" | grep '^NEW_' > "${COMPARE_DIR}/new_samples.txt" || true
bcftools query -l "${COMPARE_DIR}/LEPC_new_vs_old_joint.vcf.gz" | grep '^OLD_' > "${COMPARE_DIR}/old_samples.txt" || true

bcftools stats -S "${COMPARE_DIR}/new_samples.txt" "${COMPARE_DIR}/LEPC_new_vs_old_joint.vcf.gz" > "${COMPARE_DIR}/new_stats.txt"
bcftools stats -S "${COMPARE_DIR}/old_samples.txt" "${COMPARE_DIR}/LEPC_new_vs_old_joint.vcf.gz" > "${COMPARE_DIR}/old_stats.txt"

{
    echo -e "metric\tNEW\tOLD"
    # SN-section counts: label text lives in column 3, value in the last
    # column — straightforward substring match on each label.
    for label in "number of SNPs:" "number of indels:" "number of multiallelic sites:"; do
        new_val=$(grep "^SN" "${COMPARE_DIR}/new_stats.txt" | grep -F "${label}" | awk -F'\t' '{print $NF}' | head -n1)
        old_val=$(grep "^SN" "${COMPARE_DIR}/old_stats.txt" | grep -F "${label}" | awk -F'\t' '{print $NF}' | head -n1)
        echo -e "${label}\t${new_val:-NA}\t${old_val:-NA}"
    done
    # ts/tv ratio: a separate TSTV data row (not an SN line), format
    # TSTV, id, ts, tv, ts/tv, ts(1st ALT), tv(1st ALT), ts/tv(1st ALT) —
    # ratio is column 5, not the last column.
    new_tstv=$(grep "^TSTV" "${COMPARE_DIR}/new_stats.txt" | awk -F'\t' '{print $5}' | head -n1)
    old_tstv=$(grep "^TSTV" "${COMPARE_DIR}/old_stats.txt" | awk -F'\t' '{print $5}' | head -n1)
    echo -e "ts/tv\t${new_tstv:-NA}\t${old_tstv:-NA}"
} > "${COMPARE_DIR}/new_vs_old_summary.tsv"

echo "  Per-batch stats : ${COMPARE_DIR}/{new,old}_stats.txt (full bcftools stats output)"
echo "  Quick summary   : ${COMPARE_DIR}/new_vs_old_summary.tsv"
echo "  Full QC (per-sample depth, duplication, etc. for both batches side by side):"
echo "    ${OUTDIR}/multiqc/multiqc_report.html"
echo ""
echo ">>> End time: $(date)"
