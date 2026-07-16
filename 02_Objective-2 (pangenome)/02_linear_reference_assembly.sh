#!/bin/bash
# =============================================================================
# SLURM ARRAY JOB: DE NOVO GENOME ASSEMBLY — hifiasm + yahs
# Step 03 — requires 02_download_chicken_ref.sh to have completed first.
# n=23 samples across 3 species (Sharp-tailed Grouse, Lesser Prairie-Chicken,
# Greater Prairie-Chicken). One array task per sample.
#
# Per sample:
#   1. hifiasm assembles PacBio HiFi reads, phased into two haplotypes using
#      the sample's Hi-C reads (--h1/--h2 Hi-C integration mode). The 3
#      samples that also have Oxford Nanopore ultra-long reads pass them in
#      via --ul for extra contiguity.
#   2. Each haplotype's contigs are Hi-C scaffolded with yahs: Hi-C reads are
#      re-aligned to that haplotype's contigs (bwa mem, Arima-style -5SP),
#      filtered to MAPQ >= 20 for confident placement, then yahs scaffolds
#      using those alignments.
#   3. RagTag orders and orients the yahs scaffolds into pseudo-chromosomes
#      against the chicken (Gallus gallus, GRCg7b) reference.
#   4. Liftoff transfers chicken gene annotations onto each pseudo-chromosome
#      assembly.
#
# INPUT MANIFEST (tab-separated, header required, see assembly_manifest.tsv):
#   sample_id  species  hifi_fastqs  hic_r1  hic_r2  ont_ul
#     - hifi_fastqs : one or more HiFi fastq.gz paths, comma-separated
#                     (multiple SMRT cells are common)
#     - hic_r1/hic_r2 : paired-end Hi-C fastq.gz
#     - ont_ul      : ONT ultra-long fastq.gz, or "NA" for the 20 samples
#                     that don't have it
#
# USAGE:
#   Run 02_download_chicken_ref.sh first — it fetches the reference this
#   script needs and prints the exact sbatch command to submit this array.
#   Or submit manually:
#
#     N=$(grep -vc '^#' assembly_manifest.tsv)   # includes header, so N = samples+1
#     N=$((N - 1))
#     sbatch --array=0-$((N-1))%6 03_genome_assembly_array.sh
#   (the %6 caps concurrent tasks — hifiasm is memory-hungry; tune to your
#   cluster's fair-share policy and node memory)
# =============================================================================
#SBATCH --job-name=grouse_asm
#SBATCH --output=logs/%x_%A_%a.out
#SBATCH --error=logs/%x_%A_%a.err
#SBATCH -A fnrdewoody
#SBATCH -t 10-00:00:00
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=48
#SBATCH --mem=250G
#SBATCH -p cpu
#SBATCH --mail-type=BEGIN,END,FAIL
#SBATCH --mail-user=${USER}@purdue.edu

# =============================================================================
# ENVIRONMENT SETUP
# =============================================================================
set -euo pipefail

ml biocontainers
ml hifiasm
ml bwa
ml samtools/1.22.1
ml yahs # module availability varies by cluster — may need to build from source
         # (https://github.com/c-zhou/yahs) if `ml yahs` fails here
ml ragtag   # or `pip install ragtag` / conda if no module exists
ml liftoff  # or `pip install liftoff` / conda if no module exists

# =============================================================================
# USER-DEFINED VARIABLES
# =============================================================================
MANIFEST="${SLURM_SUBMIT_DIR}/assembly_manifest.tsv"

PROJECT_DIR="${CLUSTER_SCRATCH}/GROUSE_ASM"
REF_DIR="${PROJECT_DIR}/ref"
ASM_DIR="${PROJECT_DIR}/hifiasm"
SCAFFOLD_DIR="${PROJECT_DIR}/yahs"
RAGTAG_DIR="${PROJECT_DIR}/ragtag"
LIFTOFF_DIR="${PROJECT_DIR}/liftoff"
FINAL_DIR="${PROJECT_DIR}/final"

# Chicken (Gallus gallus) RefSeq reference — downloaded by
# 02_download_chicken_ref.sh (GCF_016699485.2 / GRCg7b); used here by RagTag
# to orient/order scaffolds into pseudo-chromosomes and by Liftoff as the
# annotation source.
CHICKEN_ASM_DIR="GCF_016699485.2_bGalGal1.mat.broiler.GRCg7b"
CHICKEN_FASTA="${REF_DIR}/${CHICKEN_ASM_DIR}_genomic.fna"
CHICKEN_GFF="${REF_DIR}/${CHICKEN_ASM_DIR}_genomic.gff"

# Minimum mapping quality for Hi-C alignments passed to yahs — reads below
# this are dropped so scaffolding only relies on confidently-placed pairs.
HIC_MIN_MAPQ=20

THREADS=$SLURM_CPUS_PER_TASK

mkdir -p logs "$ASM_DIR" "$SCAFFOLD_DIR" "$RAGTAG_DIR" "$LIFTOFF_DIR" "$FINAL_DIR" "$REF_DIR"

# =============================================================================
# VALIDATE CHICKEN REFERENCE
# Downloaded by 02_download_chicken_ref.sh — must be run before this script.
# =============================================================================
if [[ ! -f "$CHICKEN_FASTA" || ! -f "$CHICKEN_GFF" ]]; then
    echo "ERROR: Chicken reference not found at ${REF_DIR}"
    echo "Run 02_download_chicken_ref.sh first, then submit this array."
    exit 1
fi

# =============================================================================
# RESOLVE SAMPLE FOR THIS ARRAY TASK
# =============================================================================
if [[ -z "${SLURM_ARRAY_TASK_ID:-}" ]]; then
    echo "ERROR: SLURM_ARRAY_TASK_ID is not set."
    echo "Submit with: sbatch --array=0-N 03_genome_assembly_array.sh"
    exit 1
fi

if [[ ! -f "$MANIFEST" ]]; then
    echo "ERROR: Manifest not found: ${MANIFEST}"
    exit 1
fi

mapfile -t ROWS < <(grep -v '^#' "$MANIFEST" | tail -n +2)
LINE="${ROWS[$SLURM_ARRAY_TASK_ID]:-}"

if [[ -z "$LINE" ]]; then
    echo "ERROR: No manifest row at index ${SLURM_ARRAY_TASK_ID}"
    exit 1
fi

IFS=$'\t' read -r SAMPLE SPECIES HIFI_READS HIC_R1 HIC_R2 ONT_UL <<< "$LINE"

if [[ -z "$SAMPLE" || -z "$HIFI_READS" || -z "$HIC_R1" || -z "$HIC_R2" ]]; then
    echo "ERROR: Malformed manifest row at index ${SLURM_ARRAY_TASK_ID}: ${LINE}"
    exit 1
fi

echo ">>> Array task ${SLURM_ARRAY_TASK_ID} -> sample: ${SAMPLE} (${SPECIES})"
echo ">>> HiFi reads : ${HIFI_READS}"
echo ">>> Hi-C reads : ${HIC_R1} / ${HIC_R2}"
echo ">>> ONT UL     : ${ONT_UL}"
echo ">>> Running on : $(hostname)"
echo ">>> CPUs       : ${THREADS}"
echo ">>> Start time : $(date)"

IFS=',' read -ra HIFI_ARR <<< "$HIFI_READS"
for F in "${HIFI_ARR[@]}" "$HIC_R1" "$HIC_R2"; do
    if [[ ! -f "$F" ]]; then
        echo "ERROR: Input file not found: ${F}"
        exit 1
    fi
done

HAS_UL=false
if [[ -n "$ONT_UL" && "$ONT_UL" != "NA" ]]; then
    if [[ ! -f "$ONT_UL" ]]; then
        echo "ERROR: ONT UL file not found: ${ONT_UL}"
        exit 1
    fi
    HAS_UL=true
fi

# =============================================================================
# STEP 1: hifiasm assembly (Hi-C-phased, +UL when available)
# Produces haplotype-resolved primary contigs:
#   ${OUT_PREFIX}.hic.hap1.p_ctg.gfa
#   ${OUT_PREFIX}.hic.hap2.p_ctg.gfa
# =============================================================================
echo ">>> Step 1: hifiasm assembly"

OUT_PREFIX="${ASM_DIR}/${SAMPLE}"

HIFIASM_CMD=(hifiasm -o "$OUT_PREFIX" -t "$THREADS" --h1 "$HIC_R1" --h2 "$HIC_R2")
if [[ "$HAS_UL" == true ]]; then
    echo "  Ultra-long ONT reads detected — adding --ul"
    HIFIASM_CMD+=(--ul "$ONT_UL")
fi
HIFIASM_CMD+=("${HIFI_ARR[@]}")

echo "  ${HIFIASM_CMD[*]}"
"${HIFIASM_CMD[@]}"

# =============================================================================
# STEP 2: Convert phased contig GFAs to FASTA, then Hi-C scaffold each
# haplotype with yahs (Arima-style bwa mem -5SP mapping -> name-sorted BAM).
# =============================================================================
for HAP in hap1 hap2; do
    echo ">>> Step 2 (${HAP}): GFA -> FASTA"

    GFA="${OUT_PREFIX}.hic.${HAP}.p_ctg.gfa"
    CONTIGS="${ASM_DIR}/${SAMPLE}.${HAP}.contigs.fa"

    if [[ ! -f "$GFA" ]]; then
        echo "ERROR: Expected hifiasm output not found: ${GFA}"
        exit 1
    fi

    awk '/^S/{print ">"$2"\n"$3}' "$GFA" > "$CONTIGS"
    samtools faidx "$CONTIGS"

    echo ">>> Step 3 (${HAP}): Align Hi-C reads to ${HAP} contigs (MAPQ >= ${HIC_MIN_MAPQ})"

    bwa index "$CONTIGS"

    HIC_BAM="${SCAFFOLD_DIR}/${SAMPLE}.${HAP}.hic2contigs.bam"
    bwa mem -5SP -t "$THREADS" "$CONTIGS" "$HIC_R1" "$HIC_R2" \
        | samtools view -@ "$THREADS" -buS -q "$HIC_MIN_MAPQ" - \
        | samtools sort -@ "$THREADS" -n -o "$HIC_BAM" -

    echo ">>> Step 4 (${HAP}): Hi-C scaffolding with yahs"

    YAHS_PREFIX="${SCAFFOLD_DIR}/${SAMPLE}.${HAP}"
    yahs "$CONTIGS" "$HIC_BAM" -o "$YAHS_PREFIX"

    SCAFFOLDS="${YAHS_PREFIX}_scaffolds_final.fa"
    if [[ ! -f "$SCAFFOLDS" ]]; then
        echo "ERROR: yahs did not produce expected output: ${SCAFFOLDS}"
        exit 1
    fi

    echo ">>> Step 5 (${HAP}): RagTag pseudo-chromosome scaffolding vs chicken reference"

    RAGTAG_OUT="${RAGTAG_DIR}/${SAMPLE}.${HAP}"
    ragtag.py scaffold \
        -o "$RAGTAG_OUT" \
        -t "$THREADS" \
        -u \
        "$CHICKEN_FASTA" \
        "$SCAFFOLDS"

    PSEUDO_CHR="${RAGTAG_OUT}/ragtag.scaffold.fasta"
    if [[ ! -f "$PSEUDO_CHR" ]]; then
        echo "ERROR: RagTag did not produce expected output: ${PSEUDO_CHR}"
        exit 1
    fi

    FINAL_FASTA="${FINAL_DIR}/${SPECIES}_${SAMPLE}_${HAP}.pseudo_chr.fasta"
    cp "$PSEUDO_CHR" "$FINAL_FASTA"
    samtools faidx "$FINAL_FASTA"
    echo "  Final ${HAP} pseudo-chromosome assembly: ${FINAL_FASTA}"

    echo ">>> Step 6 (${HAP}): Liftoff chicken annotation onto pseudo-chromosomes"

    LIFTOFF_GFF="${FINAL_DIR}/${SPECIES}_${SAMPLE}_${HAP}.liftoff.gff3"
    LIFTOFF_UNMAPPED="${LIFTOFF_DIR}/${SAMPLE}.${HAP}.unmapped_features.txt"
    LIFTOFF_INTERMEDIATE="${LIFTOFF_DIR}/${SAMPLE}.${HAP}_intermediate"

    mkdir -p "$LIFTOFF_INTERMEDIATE"

    liftoff \
        -g "$CHICKEN_GFF" \
        -o "$LIFTOFF_GFF" \
        -u "$LIFTOFF_UNMAPPED" \
        -dir "$LIFTOFF_INTERMEDIATE" \
        -p "$THREADS" \
        "$FINAL_FASTA" \
        "$CHICKEN_FASTA"

    echo "  Liftoff annotation : ${LIFTOFF_GFF}"
    echo "  Unmapped features  : ${LIFTOFF_UNMAPPED}"
done

echo ">>> Sample ${SAMPLE} (${SPECIES}) complete — $(date)"
