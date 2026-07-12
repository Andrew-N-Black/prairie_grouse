#!/bin/bash
# =============================================================================
# SLURM JOB SUBMISSION: nf-core/sarek — LEPC GERMLINE VARIANT CALLING
# ~100 raw paired-end FASTQ samples -> GATK4 alignment/BQSR -> joint HaplotypeCaller
# genotyping -> SnpEff annotation -> final joint VCF
#
# This job only orchestrates Nextflow; the heavy per-sample work is farmed out
# as individual SLURM jobs via the 'slurm' executor profile defined in
# nextflow.config, so the head job itself needs modest resources.
# =============================================================================
#SBATCH --job-name=lepc_sarek
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
#SBATCH --mail-user=blackan@purdue.edu

# =============================================================================
# ENVIRONMENT SETUP
# =============================================================================
set -euo pipefail
mkdir -p logs

# =============================================================================
# USER-DEFINED VARIABLES
# =============================================================================
PROJECT="LEPC"

# Reference: Tympanuchus pallidicinctus (Lesser Prairie-Chicken), RefSeq
# GCF_026119805.1 / assembly pur_lepc_1.0 (Behrens et al. 2023, GBE)
ASSEMBLY_ACC="GCF_026119805.1"
ASSEMBLY_NAME="pur_lepc_1.0"
ASSEMBLY_DIR="${ASSEMBLY_ACC}_${ASSEMBLY_NAME}"
NCBI_FTP_BASE="https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/026/119/805/${ASSEMBLY_DIR}"
REFERENCE="${ASSEMBLY_DIR}_genomic.fa.gz"
ANNOTATION="${ASSEMBLY_DIR}_genomic.gff.gz"
GENOME_VERSION="Tympanuchus_pallidicinctus_pur_lepc_1.0"

# Paths derived from project root
PROJECT_DIR="${CLUSTER_SCRATCH}/${PROJECT}"
REF_DIR="${PROJECT_DIR}/ref"
RAW_DIR="${PROJECT_DIR}/raw"
OUTDIR="${PROJECT_DIR}/nf-out"
SNPEFF_CACHE="${REF_DIR}/snpeff_cache"

REF_FASTA="${REF_DIR}/${REFERENCE}"
SAMPLESHEET="${PROJECT_DIR}/samplesheet.csv"

# nf-core/sarek settings
SAREK_VERSION="3.9.0"

echo ">>> 01_obj.sh — LEPC nf-core/sarek germline pipeline"
echo ">>> Project   : ${PROJECT}"
echo ">>> Reference : ${ASSEMBLY_ACC} (${ASSEMBLY_NAME})"
echo ">>> Output dir: ${OUTDIR}"
echo ">>> Node      : $(hostname)"
echo ">>> Start time: $(date)"

mkdir -p "${REF_DIR}" "${RAW_DIR}" "${OUTDIR}" "${SNPEFF_CACHE}"

# =============================================================================
# STEP 1: Download reference genome + GFF3 annotation (skipped if present)
# =============================================================================
echo ">>> Step 1: Reference genome + annotation"

cd "${REF_DIR}"

if [[ ! -f "${REFERENCE}" ]]; then
    echo "  ${REFERENCE} not found — downloading from NCBI"
    wget -q "${NCBI_FTP_BASE}/${ASSEMBLY_DIR}_genomic.fna.gz" -O "${REFERENCE}"
    echo "  Download complete"
else
    echo "  ${REFERENCE} already present — skipping download"
fi

if [[ ! -f "${ANNOTATION}" ]]; then
    echo "  ${ANNOTATION} not found — downloading from NCBI"
    wget -q "${NCBI_FTP_BASE}/${ASSEMBLY_DIR}_genomic.gff.gz" -O "${ANNOTATION}"
    echo "  Download complete"
else
    echo "  ${ANNOTATION} already present — skipping download"
fi

# =============================================================================
# STEP 2: Build a custom SnpEff database from the RefSeq GFF3
# LEPC has no prebuilt SnpEff/VEP cache, so one is built locally from the
# reference + annotation and passed into sarek's annotation step.
# -noCheckCds/-noCheckProtein tolerate the minor CDS/protein inconsistencies
# common in draft RefSeq annotations for non-model species.
# =============================================================================
echo ">>> Step 2: Build custom SnpEff database (${GENOME_VERSION})"

ml biocontainers
ml snpeff

SNPEFF_DB_DIR="${SNPEFF_CACHE}/${GENOME_VERSION}"

if [[ ! -f "${SNPEFF_DB_DIR}/snpEffectPredictor.bin" ]]; then
    mkdir -p "${SNPEFF_DB_DIR}"
    cp "${REF_DIR}/${REFERENCE}" "${SNPEFF_DB_DIR}/sequences.fa.gz"
    cp "${REF_DIR}/${ANNOTATION}" "${SNPEFF_DB_DIR}/genes.gff.gz"

    cat > "${SNPEFF_CACHE}/snpEff.config" <<EOF
${GENOME_VERSION}.genome : Tympanuchus_pallidicinctus
EOF

    snpEff build -gff3 -noCheckCds -noCheckProtein \
        -c "${SNPEFF_CACHE}/snpEff.config" \
        -dataDir "${SNPEFF_CACHE}" \
        -v "${GENOME_VERSION}"

    echo "  SnpEff database built at ${SNPEFF_DB_DIR}"
else
    echo "  SnpEff database already built — skipping"
fi

# =============================================================================
# STEP 3: Auto-generate the sarek samplesheet from paired-end FASTQs in RAW_DIR
# Handles *_R1_001/_R2_001, *_R1/_R2, and *_1/_2 naming conventions.
# Each pair is treated as one germline sample (patient == sample, lane 1).
# =============================================================================
echo ">>> Step 3: Building sarek samplesheet from ${RAW_DIR}"

echo "patient,sample,lane,fastq_1,fastq_2" > "${SAMPLESHEET}"

N_PAIRS=0
while IFS= read -r R1; do
    BASENAME=$(basename "${R1}")
    case "${BASENAME}" in
        *_R1_001.fastq.gz) SAMPLE="${BASENAME%_R1_001.fastq.gz}"; R2="${R1/_R1_001.fastq.gz/_R2_001.fastq.gz}" ;;
        *_R1.fastq.gz)     SAMPLE="${BASENAME%_R1.fastq.gz}";     R2="${R1/_R1.fastq.gz/_R2.fastq.gz}" ;;
        *_1.fastq.gz)      SAMPLE="${BASENAME%_1.fastq.gz}";      R2="${R1/_1.fastq.gz/_2.fastq.gz}" ;;
        *) echo "  WARNING: unrecognized naming pattern for ${BASENAME} — skipping" >&2; continue ;;
    esac

    if [[ ! -f "${R2}" ]]; then
        echo "  WARNING: mate not found for ${R1} — skipping" >&2
        continue
    fi

    SAMPLE=$(echo "${SAMPLE}" | tr -c 'A-Za-z0-9_-' '_')
    echo "${SAMPLE},${SAMPLE},1,${R1},${R2}" >> "${SAMPLESHEET}"
    N_PAIRS=$((N_PAIRS + 1))
done < <(find "${RAW_DIR}" -maxdepth 1 -type f \( -name "*_R1_001.fastq.gz" -o -name "*_R1.fastq.gz" -o -name "*_1.fastq.gz" \) | sort)

echo "  Wrote ${N_PAIRS} sample pairs to ${SAMPLESHEET}"

if [[ "${N_PAIRS}" -eq 0 ]]; then
    echo "ERROR: no paired-end FASTQs found in ${RAW_DIR}" >&2
    exit 1
fi

# =============================================================================
# STEP 4: nf-core/sarek
# fastp adapter clipping -> BWA-MEM2 alignment -> GATK4 MarkDuplicates + BQSR
# -> GATK4 HaplotypeCaller (per-sample GVCF) -> joint genotyping across all
# ~100 samples -> SnpEff functional annotation -> final joint VCF.
#
# --trim_fastq       : fastp-based adapter/quality clipping
# --tools             : GATK HaplotypeCaller + SnpEff annotation
# --joint_germline    : cohort-wide joint genotyping (GenomicsDBImport/GenotypeGVCFs)
# --save_reference    : caches the BWA/dict/fai indices sarek builds from --fasta
# =============================================================================
echo ">>> Step 4: nf-core/sarek v${SAREK_VERSION}"

module load nf-core/2.11.1

cd "${PROJECT_DIR}"

nextflow run nf-core/sarek \
    -r "${SAREK_VERSION}" \
    -c "${PROJECT_DIR}/nextflow.config" \
    -profile singularity,slurm \
    --input "${SAMPLESHEET}" \
    --outdir "${OUTDIR}" \
    --fasta "${REF_FASTA}" \
    --trim_fastq \
    --tools haplotypecaller,snpeff \
    --joint_germline \
    --snpeff_cache "${SNPEFF_CACHE}" \
    --snpeff_db "${GENOME_VERSION}" \
    --save_reference \
    --email "${USER}@purdue.edu"

# =============================================================================
# STEP 5: Stage the final joint, annotated VCF to a predictable location
# sarek's exact output path can shift slightly between versions, so this
# searches the results dir rather than hardcoding one; verify against
# ${OUTDIR}/pipeline_info/execution_report.html if nothing is found.
# =============================================================================
echo ">>> Step 5: Locating final joint VCF"

mkdir -p "${PROJECT_DIR}/final"

FINAL_VCF=$(find "${OUTDIR}/annotation" -iname "*joint_germline*ann.vcf.gz" 2>/dev/null | head -n1)
if [[ -z "${FINAL_VCF}" ]]; then
    FINAL_VCF=$(find "${OUTDIR}/variant_calling" -ipath "*joint_variant_calling*" -name "*.vcf.gz" 2>/dev/null | head -n1)
fi

if [[ -n "${FINAL_VCF}" ]]; then
    cp "${FINAL_VCF}" "${PROJECT_DIR}/final/LEPC_joint_annotated.vcf.gz"
    echo "  Final VCF: ${PROJECT_DIR}/final/LEPC_joint_annotated.vcf.gz"
else
    echo "  WARNING: could not auto-locate the final VCF — check ${OUTDIR}/variant_calling and ${OUTDIR}/annotation manually"
fi

echo ""
echo ">>> End time: $(date)"
