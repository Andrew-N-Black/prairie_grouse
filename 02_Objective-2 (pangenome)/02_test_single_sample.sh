#!/bin/bash
# =============================================================================
# SLURM JOB: SINGLE-SAMPLE TEST OF 03_genome_assembly_array.sh
# Runs the exact same per-sample pipeline as the array script (hifiasm ->
# yahs -> Juicebox export -> RagTag -> Liftoff -> QC), but for ONE sample
# defined via environment variables instead of reading assembly_manifest.tsv
# by array index. Use this to confirm the pipeline actually works on a real
# sample before committing to the full array (which is expensive to debug
# once 23 tasks are all failing the same way).
#
# Requires 02_download_chicken_ref.sh to have completed first (same as the
# array script) — this reuses the same chicken reference, yahs conda env,
# and juicer_tools.jar rather than rebuilding them, since those are shared,
# expensive, sample-independent assets. Only the OUTPUT directories are
# separated (under test_single_sample/) so this run can't collide with or
# be mistaken for real array output later.
#
# USAGE:
#   export SAMPLE="your_sample_id"
#   export SPECIES="your_species_label"          # matches assembly_manifest.tsv convention
#   export HIFI_READS="/path/to/hifi_1.fastq.gz,/path/to/hifi_2.fastq.gz"
#   export HIC_R1="/path/to/hic_R1.fastq.gz"
#   export HIC_R2="/path/to/hic_R2.fastq.gz"
#   export ONT_UL="NA"                            # or a real path if this sample has UL reads
#   sbatch 03_test_single_sample.sh
#
#   (sbatch passes through the submitting shell's environment by default,
#   so the exports above reach the job. Equivalently, in one line:
#   sbatch --export=ALL,SAMPLE=...,SPECIES=...,HIFI_READS=...,HIC_R1=...,HIC_R2=...,ONT_UL=... \
#       03_test_single_sample.sh)
# =============================================================================
#SBATCH --job-name=grouse_asm_test
#SBATCH --output=logs/%x_%j.out
#SBATCH --error=logs/%x_%j.err
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
# No yahs module on this cluster — installed below via conda/bioconda into a
# dedicated env under PROJECT_DIR (see YAHS_ENV_DIR further down). That
# package also provides the `juicer` binary used in Step 5 (yahs's own JBAT
# pre-processor, bundled in the same repo — NOT Aidenlab's official juicer
# pipeline of the same name).
ml anaconda/2025.12-py313
ml java # for juicer_tools.jar (Step 5)
ml ragtag   # or `pip install ragtag` / conda if no module exists
ml liftoff  # or `pip install liftoff` / conda if no module exists
ml minimap2
ml pretextmap       # or build from https://github.com/sanger-tol/PretextMap
ml pretextsnapshot   # or build from https://github.com/sanger-tol/PretextSnapshot
ml python            # needs matplotlib for the orientation dotplot

# =============================================================================
# SAMPLE VARIABLES — set these via `export` before calling sbatch (see
# USAGE above). Fails immediately and clearly if any are missing, rather
# than proceeding with an empty/wrong value.
# =============================================================================
: "${SAMPLE:?Set SAMPLE before submitting, e.g. export SAMPLE=your_sample_id}"
: "${SPECIES:?Set SPECIES before submitting, e.g. export SPECIES=lesser_prairie-chicken}"
: "${HIFI_READS:?Set HIFI_READS before submitting (comma-separated if multiple files)}"
: "${HIC_R1:?Set HIC_R1 before submitting}"
: "${HIC_R2:?Set HIC_R2 before submitting}"
: "${ONT_UL:=NA}"   # optional — defaults to NA (no ultra-long reads) if unset

# =============================================================================
# USER-DEFINED VARIABLES
# Same shared, expensive, sample-independent assets as the array script
# (chicken reference, yahs conda env, juicer_tools.jar) — reused as-is, not
# rebuilt. Only output paths differ (nested under test_single_sample/).
# =============================================================================
PROJECT_DIR="${CLUSTER_SCRATCH}/GROUSE_ASM"
REF_DIR="${PROJECT_DIR}/ref"

TEST_DIR="${PROJECT_DIR}/test_single_sample"
ASM_DIR="${TEST_DIR}/hifiasm"
SCAFFOLD_DIR="${TEST_DIR}/yahs"
RAGTAG_DIR="${TEST_DIR}/ragtag"
LIFTOFF_DIR="${TEST_DIR}/liftoff"
FINAL_DIR="${TEST_DIR}/final"
QC_DIR="${TEST_DIR}/qc"

CHICKEN_ASM_DIR="GCF_016699485.2_bGalGal1.mat.broiler.GRCg7b"
CHICKEN_FASTA="${REF_DIR}/${CHICKEN_ASM_DIR}_genomic.fna"
CHICKEN_GFF="${REF_DIR}/${CHICKEN_ASM_DIR}_genomic.gff"

HIC_MIN_MAPQ=20

JUICER_TOOLS_JAR="${PROJECT_DIR}/juicer_tools.jar"

YAHS_VERSION="1.2.2"
YAHS_ENV_DIR="${PROJECT_DIR}/conda_envs/yahs-${YAHS_VERSION}"
YAHS_BIN="${YAHS_ENV_DIR}/bin/yahs"
JUICER_BIN="${YAHS_ENV_DIR}/bin/juicer"

THREADS=$SLURM_CPUS_PER_TASK

mkdir -p logs "$ASM_DIR" "$SCAFFOLD_DIR" "$RAGTAG_DIR" "$LIFTOFF_DIR" "$FINAL_DIR" "$QC_DIR"

# =============================================================================
# INSTALL yahs (once — skipped if the real 03/02 setup already built it,
# since YAHS_ENV_DIR points at the same shared conda env)
# =============================================================================
if [[ ! -x "$YAHS_BIN" || ! -x "$JUICER_BIN" ]]; then
    echo ">>> Installing yahs v${YAHS_VERSION} via conda (bioconda)"
    conda create --yes --prefix "$YAHS_ENV_DIR" -c bioconda -c conda-forge "yahs=${YAHS_VERSION}"
fi

if [[ ! -x "$YAHS_BIN" || ! -x "$JUICER_BIN" ]]; then
    echo "ERROR: yahs install failed — expected binaries not found in ${YAHS_ENV_DIR}/bin"
    exit 1
fi

# =============================================================================
# VALIDATE CHICKEN REFERENCE + juicer_tools.jar
# =============================================================================
if [[ ! -f "$CHICKEN_FASTA" || ! -f "$CHICKEN_GFF" ]]; then
    echo "ERROR: Chicken reference not found at ${REF_DIR}"
    echo "Run 02_download_chicken_ref.sh first, then resubmit."
    exit 1
fi

if [[ ! -f "$JUICER_TOOLS_JAR" ]]; then
    echo "ERROR: juicer_tools.jar not found at ${JUICER_TOOLS_JAR}"
    echo "Download it once from https://github.com/aidenlab/juicer/wiki/Download"
    exit 1
fi

# =============================================================================
# VALIDATE SAMPLE INPUT FILES
# =============================================================================
echo ">>> TEST RUN — single sample: ${SAMPLE} (${SPECIES})"
echo ">>> HiFi reads : ${HIFI_READS}"
echo ">>> Hi-C reads : ${HIC_R1} / ${HIC_R2}"
echo ">>> ONT UL     : ${ONT_UL}"
echo ">>> Output dir : ${TEST_DIR}"
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
# haplotype with yahs (bwa mem -5SP split-read mapping -> name-sorted BAM;
# this mapping strategy is vendor-agnostic, not specific to any Hi-C kit).
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
    "$YAHS_BIN" "$CONTIGS" "$HIC_BAM" -o "$YAHS_PREFIX"

    SCAFFOLDS="${YAHS_PREFIX}_scaffolds_final.fa"
    if [[ ! -f "$SCAFFOLDS" ]]; then
        echo "ERROR: yahs did not produce expected output: ${SCAFFOLDS}"
        exit 1
    fi

    echo ">>> Step 5 (${HAP}): Build .hic file for Juicebox visualization of raw yahs scaffolding"
    echo "    (this is yahs's own Hi-C scaffolding, BEFORE RagTag reorders it against"
    echo "     chicken — open in Juicebox to inspect/manually curate the raw Hi-C joins)"

    JBAT_PREFIX="${SCAFFOLD_DIR}/${SAMPLE}.${HAP}_JBAT"

    # yahs ships its own `juicer` binary (built alongside `yahs` from the same
    # source repo: https://github.com/c-zhou/yahs) — distinct from Aidenlab's
    # official juicer pipeline. `-a` targets Juicebox Assembly Tools (JBAT) output.
    "$JUICER_BIN" pre -a -o "$JBAT_PREFIX" \
        "${YAHS_PREFIX}.bin" "${YAHS_PREFIX}_scaffolds_final.agp" "${CONTIGS}.fai" \
        > "${JBAT_PREFIX}.log" 2>&1

    # juicer_tools pre needs a chrom-sizes-style input; yahs's own log records
    # the exact sizes to use (PRE_C_SIZE lines) — see yahs README for this idiom.
    java -Xmx32G -jar "$JUICER_TOOLS_JAR" pre \
        "${JBAT_PREFIX}.txt" "${JBAT_PREFIX}.hic.part" \
        <(grep PRE_C_SIZE "${JBAT_PREFIX}.log" | awk '{print $2" "$3}')
    mv "${JBAT_PREFIX}.hic.part" "${JBAT_PREFIX}.hic"

    echo "  Juicebox files: ${JBAT_PREFIX}.hic + ${JBAT_PREFIX}.assembly"
    echo "  Open both in Juicebox to view/curate; if total assembly length > 2Gb,"
    echo "  check ${JBAT_PREFIX}.log for a scale factor to set via Assembly > Set Scale"

    echo ">>> Step 6 (${HAP}): RagTag pseudo-chromosome scaffolding vs chicken reference"

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

    echo ">>> Step 7 (${HAP}): Liftoff chicken annotation onto pseudo-chromosomes"

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

    echo ">>> Step 8 (${HAP}): Validate scaffold order/orientation vs chicken (dotplot)"
    echo "    Independent check via minimap2 — does not reuse RagTag's own alignment."

    PAF="${QC_DIR}/${SAMPLE}.${HAP}.vs_chicken.paf"
    minimap2 -x asm20 -t "$THREADS" "$CHICKEN_FASTA" "$FINAL_FASTA" > "$PAF"

    DOTPLOT="${QC_DIR}/${SAMPLE}.${HAP}.orientation_dotplot.png"
    python3 - "$PAF" "$DOTPLOT" "${SPECIES} ${SAMPLE} ${HAP} vs chicken (GRCg7b)" << 'PYEOF'
import math, sys
from collections import defaultdict
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

paf_path, out_path, title = sys.argv[1:4]

# Group alignment blocks by query (assembly) sequence, and track total
# matched bases per (query, target) pair so each query can be paired with
# its best-matching chicken chromosome for plotting.
blocks_by_query = defaultdict(list)
matched_bases = defaultdict(lambda: defaultdict(int))

with open(paf_path) as fh:
    for line in fh:
        f = line.rstrip('\n').split('\t')
        if len(f) < 12:
            continue
        qname, qstart, qend, strand = f[0], int(f[2]), int(f[3]), f[4]
        tname, tstart, tend = f[5], int(f[7]), int(f[8])
        match_len = int(f[9])
        if match_len < 1000:  # drop short/spurious alignments
            continue
        blocks_by_query[qname].append((qstart, qend, tstart, tend, strand, tname))
        matched_bases[qname][tname] += match_len

if not blocks_by_query:
    print("No alignments above length threshold; skipping dotplot.")
    sys.exit(0)

primary_target = {q: max(tb.items(), key=lambda kv: kv[1])[0] for q, tb in matched_bases.items()}
queries = sorted(blocks_by_query.keys())

ncols = min(6, len(queries))
nrows = math.ceil(len(queries) / ncols)
fig, axes = plt.subplots(nrows, ncols, figsize=(3 * ncols, 3 * nrows), squeeze=False)

for i, qname in enumerate(queries):
    ax = axes[i // ncols][i % ncols]
    tgt = primary_target[qname]
    for qstart, qend, tstart, tend, strand, tname in blocks_by_query[qname]:
        if tname != tgt:
            continue
        color = 'tab:blue' if strand == '+' else 'tab:red'
        ys = (tstart, tend) if strand == '+' else (tend, tstart)
        ax.plot([qstart, qend], ys, color=color, linewidth=1.5)
    ax.set_title(qname, fontsize=7)
    ax.set_xlabel(f"vs {tgt}", fontsize=6)
    ax.tick_params(labelsize=5)

for j in range(len(queries), nrows * ncols):
    axes[j // ncols][j % ncols].axis('off')

fig.suptitle(f"{title}\nblue = + strand (expected orientation), red = - strand (inverted)", fontsize=10)
fig.tight_layout(rect=[0, 0, 1, 0.95])
fig.savefig(out_path, dpi=150)
print(f"Wrote {out_path}")
PYEOF

    echo "  Orientation dotplot: ${DOTPLOT}"

    echo ">>> Step 9 (${HAP}): Validate Hi-C scaffolding — static contact map (Pretext)"
    echo "    Hi-C reads re-aligned to the FINAL pseudo-chromosome assembly"
    echo "    (distinct from the pre-scaffolding alignment yahs used in Step 3)."

    bwa index "$FINAL_FASTA"

    FINAL_HIC_BAM="${QC_DIR}/${SAMPLE}.${HAP}.hic2final.bam"
    bwa mem -5SP -t "$THREADS" "$FINAL_FASTA" "$HIC_R1" "$HIC_R2" \
        | samtools view -@ "$THREADS" -buS -q "$HIC_MIN_MAPQ" - \
        | samtools sort -@ "$THREADS" -o "$FINAL_HIC_BAM" -
    samtools index "$FINAL_HIC_BAM"

    PRETEXT_MAP="${QC_DIR}/${SAMPLE}.${HAP}.pretext"
    samtools view -h "$FINAL_HIC_BAM" \
        | PretextMap -o "$PRETEXT_MAP" --sortby length --sortorder descend --mapq "$HIC_MIN_MAPQ"

    PretextSnapshot --map "$PRETEXT_MAP" --sequences "=full" \
        --prefix "${SAMPLE}.${HAP}." --folder "$QC_DIR"

    echo "  Hi-C contact map    : ${QC_DIR}/${SAMPLE}.${HAP}.*.png"
    echo "  (clean scaffolding shows an unbroken diagonal per chromosome;"
    echo "   off-diagonal blocks indicate misjoins or mis-ordering)"
done

echo ""
echo ">>> TEST RUN complete for ${SAMPLE} (${SPECIES}) — $(date)"
echo ">>> If this all looks right, add your real samples to assembly_manifest.tsv"
echo ">>> and submit the full array via 03_genome_assembly_array.sh."
echo ">>> (Test output lives under ${TEST_DIR} — separate from real array output,"
echo ">>> safe to delete once you've confirmed the pipeline works.)"
