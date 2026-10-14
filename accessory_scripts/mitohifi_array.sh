#!/bin/bash
# =============================================================================
# accessory_scripts/mitohifi_array.sh
# One complete, circularised, rotated and annotated mitogenome per individual,
# assembled from the HiFi reads with MitoHiFi v3, replacing the hifiasm chr_MT.
#
# WHY: the hifiasm/RagTag chr_MT sequences are 20-50 kb, but a galliform
# mitogenome is ~16.7 kb. hifiasm does not circularise, so a circular molecule
# comes out as a linear contig carrying 1.2-3 tandem copies (a concatemer), and
# long reads that carry NUMTs (nuclear mitochondrial insertions) can be pulled
# into it. Because mtDNA is one maternal haplotype, hifiasm puts it in only one
# haplotype per bird, which is why only 18 of 46 haplotypes have a chr_MT.
#
# HOW: reads are first pre-filtered to those that map to a reference mitogenome.
# MitoHiFi then:
#   - drops reads longer than --max-read-len x the reference length; NUMT-bearing
#     reads are long nuclear reads with only part of their length mitochondrial,
#     so this removes them;
#   - assembles with hifiasm and picks the contig that best matches the reference;
#   - removes the circular overlap, giving a single copy;
#   - rotates the result to start at tRNA-Phe and annotates it (MitoFinder,
#     vertebrate mitochondrial code 2).
#
# ONE-TIME SETUP (login node), from the prairie_grouse repo root:
#   bash accessory_scripts/mitohifi_array.sh setup
#
# RUN (array over manifest rows; 23 birds):
#   N=$(grep -v '^#' assembly_manifest.tsv | tail -n +2 | grep -c .)
#   sbatch --array=0-$((N-1))%8 accessory_scripts/mitohifi_array.sh
#
# SUMMARY (login node, after the array):
#   bash accessory_scripts/mitohifi_array.sh summary
# =============================================================================
#SBATCH -A dewoody
#SBATCH -p cpu
#SBATCH -J mitohifi
#SBATCH -N 1
#SBATCH -n 1
#SBATCH -c 32
#SBATCH --mem=64G
#SBATCH -t 1-00:00:00
#SBATCH -o logs/%x_%A_%a.out
#SBATCH -e logs/%x_%A_%a.err

set -euo pipefail
unset LD_PRELOAD

PROJECT_DIR="${CLUSTER_SCRATCH}/GROUSE/grouse_asm"
MANIFEST="${SLURM_SUBMIT_DIR:-$PWD}/assembly_manifest.tsv"
FILT_DIR="${PROJECT_DIR}/hifi_filtered"            # <SAMPLE>/*.filt.fastq.gz (step 02)
MITO_DIR="${PROJECT_DIR}/mitogenomes"
REF_DIR="${MITO_DIR}/reference"
SIF="${PROJECT_DIR}/containers/mitohifi.sif"
THREADS="${SLURM_CPUS_PER_TASK:-8}"
EMAIL="blackan@purdue.edu"                          # NCBI Entrez needs an address

# Reference mitogenome. A Tympanuchus genome is used if NCBI has one, otherwise
# rock ptarmigan (same subfamily, Tetraoninae).
REF_SPECIES_PRIMARY="Tympanuchus cupido"
REF_SPECIES_FALLBACK="Lagopus muta"

# Reads longer than this multiple of the reference length are dropped
# (anti-NUMT filter). 1.0 is MitoHiFi's default; relax to 1.5 only if too few
# reads remain for an assembly.
MAX_READ_LEN="${MAX_READ_LEN:-1.0}"

APPTAINER_BIN="$(command -v apptainer || command -v singularity)"
cx() { "$APPTAINER_BIN" exec --cleanenv --bind "${CLUSTER_SCRATCH}" "$SIF" "$@"; }

mkdir -p logs "$MITO_DIR" "$REF_DIR"

# ---------------------------------------------------------------- setup -------
if [[ "${1:-}" == "setup" ]]; then
    mkdir -p "$(dirname "$SIF")"
    [[ -s "$SIF" ]] || "$APPTAINER_BIN" pull "$SIF" docker://ghcr.io/marcelauliano/mitohifi:master
    cd "$REF_DIR"
    for SP in "$REF_SPECIES_PRIMARY" "$REF_SPECIES_FALLBACK"; do
        echo ">>> Looking for a complete mitogenome of ${SP}"
        if cx findMitoReference.py --species "$SP" --email "$EMAIL" \
               --outfolder "$REF_DIR" --min_length 15000 && ls "$REF_DIR"/*.fasta >/dev/null 2>&1; then
            break
        fi
    done
    REF_FA=$(ls "$REF_DIR"/*.fasta | head -n1)
    REF_GB=$(ls "$REF_DIR"/*.gb | head -n1)
    ln -sfn "$REF_FA" "$REF_DIR/ref.fasta"; ln -sfn "$REF_GB" "$REF_DIR/ref.gb"
    echo ">>> Reference: $(grep '>' "$REF_FA")"
    echo ">>> Length   : $(grep -v '>' "$REF_FA" | tr -d '\n' | wc -c) bp"
    exit 0
fi

# ---------------------------------------------------------------- summary -----
if [[ "${1:-}" == "summary" ]]; then
    printf "sample\tlength_bp\tn_genes_annotated\tcontigs_considered\n"
    for D in "$MITO_DIR"/*/; do
        S=$(basename "$D"); [[ "$S" == reference ]] && continue
        F="${D}final_mitogenome.fasta"
        if [[ -s "$F" ]]; then
            L=$(grep -v '>' "$F" | tr -d '\n' | wc -c)
            G=$( (grep -c 'gene' "${D}final_mitogenome.gff" 2>/dev/null) || echo NA)
            C=$( (tail -n +2 "${D}contigs_stats.tsv" 2>/dev/null | wc -l) || echo NA)
            printf "%s\t%s\t%s\t%s\n" "$S" "$L" "$G" "$C"
        else
            printf "%s\tFAILED\t-\t-\n" "$S"
        fi
    done
    exit 0
fi

# ---------------------------------------------------------------- per sample --
[[ -s "$SIF" && -s "$REF_DIR/ref.fasta" ]] || { echo "ERROR: run '$0 setup' first"; exit 1; }
LINE=$(grep -v '^#' "$MANIFEST" | tail -n +2 | sed -n "$((SLURM_ARRAY_TASK_ID + 1))p")
IFS=$'\t' read -r SAMPLE SPECIES _ <<< "$LINE"
SPECIES="${SPECIES//[[:space:]]/}"
OUT="${MITO_DIR}/${SPECIES}_${SAMPLE}"
mkdir -p "$OUT"; cd "$OUT"
echo ">>> ${SPECIES} ${SAMPLE} -> ${OUT}  (${THREADS} threads, max-read-len ${MAX_READ_LEN}x)"

if [[ -s final_mitogenome.fasta ]]; then echo "already done — skipping"; exit 0; fi

READS=( "${FILT_DIR}/${SAMPLE}"/*.filt.fastq.gz )
[[ -e "${READS[0]}" ]] || { echo "ERROR: no filtered HiFi reads in ${FILT_DIR}/${SAMPLE}"; exit 1; }

# 1. Keep only reads that map to the reference mitogenome (primary alignments).
#    This shrinks ~30 Gb of reads to a few thousand. NUMT reads still pass here;
#    MitoHiFi's length filter removes them in step 2.
if [[ ! -s mito_candidate_reads.fa ]]; then
    cx bash -c "minimap2 -t ${THREADS} -ax map-hifi --secondary=no '${REF_DIR}/ref.fasta' ${READS[*]} \
        | samtools fasta -F 0x904 - > mito_candidate_reads.fa"
fi
echo ">>> candidate reads: $(grep -c '>' mito_candidate_reads.fa)"

# 2. MitoHiFi (-o 2 = vertebrate mitochondrial code; -a animal)
cx mitohifi.py -r mito_candidate_reads.fa \
    -f "${REF_DIR}/ref.fasta" -g "${REF_DIR}/ref.gb" \
    -t "$THREADS" -o 2 -a animal --max-read-len "$MAX_READ_LEN"

echo ">>> final length: $(grep -v '>' final_mitogenome.fasta | tr -d '\n' | wc -c) bp"
