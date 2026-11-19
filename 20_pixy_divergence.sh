#!/bin/bash
# =============================================================================
# SLURM JOB SUBMISSION: DXY, FST AND PI ACROSS THE GENOME WITH PIXY
# Step 07 — requires step 03 (pangenome/chr_*/*.final.og) and step 05
# (structural_variants/chr_*/ref.txt: the reference path of every chromosome).
#
# pixy needs an all-sites VCF (invariant sites included) in one linear coordinate
# system, which a graph cannot provide directly. So:
#   1. ref:     the reference path of every chromosome (the same paths step 05
#               used, so windows share coordinates with the SV and Z-inversion
#               results) is extracted from the graph into reference.fa. chr_W
#               (STGR path) is included as a decoy so female W reads do not map
#               to the Z gametologs; it is not analyzed.
#   2. map:     each bird's HiFi reads -> reference.fa (minimap2 -x map-hifi), array
#   3. call:    per 5-Mb chunk (array): bcftools mpileup/call on all 23 birds,
#               all sites kept; invariant and SNP sites filtered separately (pixy
#               recommendation), genotypes outside MIN_DP..MAX_DP set missing;
#               indels and multiallelic sites dropped; then pixy pi/dxy/fst in
#               WIN windows for STGR, GRPC and LEPC
#   4. summary: 07_pixy_summary.py: per-chromosome and genome-wide values after
#               excluding low-site and divergence-outlier windows, Z vs autosomes with
#               1-Mb block-bootstrap CIs, pi per species vs the neutral Z/A range
#               (0.5625-0.75), genome scan; 07_fig7_diversity_divergence.py: Figure 7
#
# SETUP (once, on a login node)
#   module load conda
#   conda create -y -n pixy -c conda-forge -c bioconda pixy bcftools samtools htslib minimap2 pandas matplotlib
#   (the script puts ${PIXY_ENV}/bin on PATH itself; never overwrite PATH)
#
# USAGE (from the repo root)
#   bash 07_pixy_divergence.sh samples      # writes samples.tsv from HIFI_DIR: CHECK IT
#   sbatch --export=ALL,STAGE=ref 07_pixy_divergence.sh
#   bash 07_pixy_divergence.sh map          # submits one mapping task per bird
#   bash 07_pixy_divergence.sh call         # submits the chunk array + summary (afterok)
#   STAGE=summary sbatch --export=ALL 07_pixy_divergence.sh      # rerun the summary only
#
# samples.tsv: sample <tab> species <tab> reads (comma-separated .fastq.gz/.fq.gz/
#   .fasta.gz or unaligned HiFi .bam) <tab> sex (F = the birds in FEMALES)
#
# NOTES
#   - Females are hemizygous for the Z but are called as diploid; their Z calls are
#     homozygous, so each female Z counts as two identical alleles. dxy (between
#     species) is unaffected in expectation; within-species pi on the Z is biased
#     slightly downward for GRPC (5 of 6 females) and STGR (1 of 10). Z chunks are
#     therefore also run with males only (17 birds; GRPC = F5597 alone), reported
#     as chr_Z_males in the summary.
#   - Reference bias check: rerun with REF_OVERRIDE=<STGR or GRPC haplotype FASTA,
#     chromosomes named chr_*> and OUT=<new dir>; dxy should barely change.
# =============================================================================
#SBATCH --job-name=pixy07
#SBATCH -A dewoody
#SBATCH -p cpu
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --mem=32G
#SBATCH -t 12:00:00
#SBATCH --output=logs/%x_%A_%a.out
#SBATCH --error=logs/%x_%A_%a.err
set -euo pipefail
unset LD_PRELOAD

PROJECT_DIR="${CLUSTER_SCRATCH}/GROUSE/grouse_asm"
GRAPH_ROOT="${PROJECT_DIR}/pangenome"
SV_ROOT="${PROJECT_DIR}/structural_variants"
OUT="${OUT:-${PROJECT_DIR}/pixy_divergence}"
PGGB_SIF="${PGGB_SIF:-${PROJECT_DIR}/containers/pggb.sif}"
PIXY_ENV="${PIXY_ENV:-/home/blackan/.conda/envs/pixy}"
HIFI_DIR="${HIFI_DIR:-${PROJECT_DIR}/raw/HiFi/bams}"   # the reads given to hifiasm
FEMALES="${FEMALES:-F5503,F5595,F5596,F5598,F5599,F5600}"   # birds whose assemblies carry a W
REF_OVERRIDE="${REF_OVERRIDE:-}"
SAMPLES="${SAMPLES:-${OUT}/samples.tsv}"

WIN="${WIN:-50000}"              # pixy window (bp)
CHUNK="${CHUNK:-5000000}"        # calling chunk (bp); must be a multiple of WIN
MIN_MQ="${MIN_MQ:-20}"
MIN_BQ="${MIN_BQ:-20}"
MIN_QUAL="${MIN_QUAL:-30}"       # SNP site QUAL
MIN_DP="${MIN_DP:-5}"            # per-genotype depth
MAX_DP="${MAX_DP:-100}"          # per-genotype depth (~2.3x the ~44x mean)
MIN_SITES="${MIN_SITES:-5000}"   # windows with fewer callable sites are excluded from all summaries
MAX_DIV="${MAX_DIV:-0.02}"       # windows with pi or dxy above this (~5x genome mean; collapsed repeats) are excluded
BOOT="${BOOT:-1000}"

STAGE="${STAGE:-${1:-}}"
REPO="${SLURM_SUBMIT_DIR:-$(pwd)}"
THREADS="${SLURM_CPUS_PER_TASK:-8}"
export PATH="${PIXY_ENV}/bin:${PATH}"

REF="${OUT}/ref/reference.fa"
mkdir -p "$OUT"/{ref,bam,vcf,pixy,summary} logs

APPTAINER_BIN="$(command -v apptainer || command -v singularity || true)"
cx() { local sif="$1"; shift; "$APPTAINER_BIN" exec --cleanenv --bind "${CLUSTER_SCRATCH}" "$sif" "$@"; }

need() { for t in "$@"; do command -v "$t" >/dev/null || { echo "ERROR: $t not found (conda env ${PIXY_ENV}, see SETUP)"; exit 1; }; done; }

case "$STAGE" in
# -----------------------------------------------------------------------------
samples)
    [[ -n "$HIFI_DIR" ]] || { echo "ERROR: set HIFI_DIR=<folder holding the HiFi reads>"; exit 1; }
    OG=$(ls "${GRAPH_ROOT}"/chr_1/*.final.og | head -n1)
    # birds and species from the PanSN path names (SPECIES.SAMPLE#hap#chr)
    cx "$PGGB_SIF" odgi paths -i "$OG" -L | cut -d'#' -f1 | sort -u > "${OUT}/birds.txt"
    tr ',' '\n' <<< "$FEMALES" > "${OUT}/females.txt"
    printf "sample\tspecies\treads\tsex\n" > "$SAMPLES"
    while IFS=. read -r SP S; do
        # top level of HIFI_DIR only (subfolders such as reseq/ are ignored); a bird's reads are
        # <ID>_*.bam / <ID>.bam (e.g. F5540_m84221_..._.hifi_reads.bc2179.bam, F5462_MERGED.bam)
        # or <ID>_*.fastq.gz / .fq.gz / .fasta.gz; .pbi indexes are not matched
        R=$(find -L "$HIFI_DIR" -maxdepth 1 -type f \( -name "${S}_*.bam" -o -name "${S}.bam" \
                -o -name "${S}_*.fastq.gz" -o -name "${S}_*.fq.gz" -o -name "${S}_*.fasta.gz" \) \
                ! -iname "*hic*" ! -iname "*hi-c*" ! -iname "*_ont*" ! -iname "*nanopore*" ! -name "*_R[12]*" 2>/dev/null | sort | paste -sd, -)
        SEX=M; grep -qx "$S" "${OUT}/females.txt" && SEX=F
        printf "%s\t%s\t%s\t%s\n" "$S" "$SP" "${R:-MISSING}" "$SEX" >> "$SAMPLES"
    done < "${OUT}/birds.txt"
    column -t -s $'\t' "$SAMPLES" | cut -c1-220
    echo
    echo "Birds: $(($(wc -l < "$SAMPLES") - 1))  females: $(wc -l < "${OUT}/females.txt")  missing reads: $(grep -c MISSING "$SAMPLES" || true)"
    echo "CHECK ${SAMPLES} (one HiFi set per bird, no Hi-C/ONT files) before the map stage."
    ;;
# -----------------------------------------------------------------------------
ref)
    need samtools minimap2
    if [[ -n "$REF_OVERRIDE" ]]; then
        echo ">>> reference from ${REF_OVERRIDE}"
        zcat -f "$REF_OVERRIDE" > "$REF"
    else
        : > "$REF"
        for d in "${SV_ROOT}"/chr_*; do
            CHR=$(basename "$d")
            [[ -s "${d}/ref.txt" ]] || { echo "  skip ${CHR}: no ref.txt"; continue; }
            P=$(cut -f1 "${d}/ref.txt")
            OG=$(ls "${GRAPH_ROOT}/${CHR}"/*.final.og | head -n1)
            cx "$PGGB_SIF" odgi paths -i "$OG" -f \
                | awk -v r="$P" -v c="$CHR" '/^>/{k = (substr($1, 2) == r); if (k) print ">" c " " r; next} k' >> "$REF"
            echo "  ${CHR} <- ${P}"
        done
    fi
    samtools faidx "$REF"
    cut -f1,2 "${REF}.fai" | awk '{s += $2} END {printf "  reference: %d sequences, %.1f Mb\n", NR, s / 1e6}'
    minimap2 -x map-hifi -t "$THREADS" -d "${OUT}/ref/reference.mmi" "$REF"
    echo ">>> ref done $(date)"
    ;;
# -----------------------------------------------------------------------------
map)
    need samtools minimap2
    if [[ -z "${SLURM_ARRAY_TASK_ID:-}" ]]; then
        [[ -s "${OUT}/ref/reference.mmi" ]] || { echo "ERROR: run the ref stage first"; exit 1; }
        grep -q MISSING "$SAMPLES" && { echo "ERROR: samples.tsv has MISSING reads"; exit 1; }
        N=$(($(wc -l < "$SAMPLES") - 1))
        sbatch --export=ALL,STAGE=map --array=1-${N}%8 -c 32 --mem=96G -t 24:00:00 -J pixy07_map "$0"
        exit 0
    fi
    IFS=$'\t' read -r S SP READS SEX < <(sed -n "$((SLURM_ARRAY_TASK_ID + 1))p" "$SAMPLES")
    BAM="${OUT}/bam/${S}.bam"
    [[ -s "${BAM}.bai" ]] && { echo "${S}: done already"; exit 0; }
    echo ">>> map ${S} (${SP}, ${SEX})  $(date)"
    {
        IFS=',' read -ra F <<< "$READS"
        for f in "${F[@]}"; do
            if [[ "$f" == *.bam ]]; then samtools fastq -@ 4 "$f" 2>/dev/null; else zcat -f "$f"; fi
        done
    } | minimap2 -ax map-hifi -t "$THREADS" --secondary=no -R "@RG\tID:${S}\tSM:${S}" "${OUT}/ref/reference.mmi" - \
      | samtools sort -@ 8 -m 3G -T "${OUT}/bam/${S}.tmp" -o "${BAM}.part" -
    mv "${BAM}.part" "$BAM"
    samtools index -@ 8 "$BAM"
    samtools coverage "$BAM" > "${OUT}/bam/${S}.coverage.tsv"
    awk 'NR > 1 && $1 != "chr_Z" && $1 != "chr_W" {b += $3 - $2 + 1; d += $7 * ($3 - $2 + 1)}
         END {printf "  mean autosomal depth: %.1fx\n", d / b}' "${OUT}/bam/${S}.coverage.tsv"
    echo ">>> done $(date)"
    ;;
# -----------------------------------------------------------------------------
call)
    need bcftools pixy tabix
    if [[ -z "${SLURM_ARRAY_TASK_ID:-}" ]]; then
        (( CHUNK % WIN == 0 )) || { echo "ERROR: CHUNK must be a multiple of WIN"; exit 1; }
        awk 'NR > 1' "$SAMPLES" | while IFS=$'\t' read -r S SP R X; do
            [[ -s "${OUT}/bam/${S}.bam.bai" ]] || { echo "ERROR: no BAM for ${S}"; exit 1; }
            echo "${OUT}/bam/${S}.bam"; done > "${OUT}/bams.txt"
        awk 'NR > 1 {print $1 "\t" $2}' "$SAMPLES" > "${OUT}/populations.txt"
        awk 'NR > 1 && $4 == "M" {print $1 "\t" $2}' "$SAMPLES" > "${OUT}/populations_males.txt"
        # chunks of every analyzed chromosome (chr_W is a decoy only)
        awk -v c="$CHUNK" '$1 != "chr_W" {for (s = 1; s <= $2; s += c) {e = s + c - 1; if (e > $2) e = $2; print $1 "\t" s "\t" e}}' \
            "${REF}.fai" > "${OUT}/chunks.txt"
        N=$(wc -l < "${OUT}/chunks.txt")
        JID=$(sbatch --parsable --export=ALL,STAGE=call --array=1-${N}%60 -c 4 --mem=16G -t 24:00:00 -J pixy07_call "$0")
        SID=$(sbatch --parsable --export=ALL,STAGE=summary --dependency=afterok:${JID} -c 4 --mem=32G -t 02:00:00 -J pixy07_sum "$0")
        echo "Submitted ${N} chunk tasks (job ${JID}); summary job ${SID} runs when all succeed"
        exit 0
    fi
    IFS=$'\t' read -r C S E < <(sed -n "${SLURM_ARRAY_TASK_ID}p" "${OUT}/chunks.txt")
    TAG="${C}_${S}"
    if [[ -s "${OUT}/pixy/${TAG}_dxy.txt" && ( "$C" != "chr_Z" || -s "${OUT}/pixy/${TAG}_males_dxy.txt" ) ]]; then
        echo "${TAG}: done already"; exit 0; fi
    T="${OUT}/vcf/${TAG}"
    echo ">>> ${C}:${S}-${E}  $(date)"
    bcftools mpileup -f "$REF" -r "${C}:${S}-${E}" -b "${OUT}/bams.txt" -X pacbio-ccs -q "$MIN_MQ" -Q "$MIN_BQ" \
        -a FORMAT/DP,FORMAT/AD --max-depth 2000 -Ou \
        | bcftools call -m -Ob -o "${T}.raw.bcf"
    # invariant sites and biallelic SNPs filtered separately; low/high-depth genotypes -> missing
    bcftools view -i 'ALT="."' "${T}.raw.bcf" -Ou \
        | bcftools +setGT -Ob -o "${T}.inv.bcf" - -- -t q -n . -i "FMT/DP<${MIN_DP} | FMT/DP>${MAX_DP}"
    bcftools view -m2 -M2 -v snps -i "QUAL>=${MIN_QUAL}" "${T}.raw.bcf" -Ou \
        | bcftools +setGT -Ob -o "${T}.snp.bcf" - -- -t q -n . -i "FMT/DP<${MIN_DP} | FMT/DP>${MAX_DP}"
    bcftools index -f "${T}.inv.bcf"; bcftools index -f "${T}.snp.bcf"
    bcftools concat -a -Oz -o "${T}.vcf.gz" "${T}.inv.bcf" "${T}.snp.bcf"
    tabix -f -p vcf "${T}.vcf.gz"
    echo "  sites: $(bcftools view -H "${T}.inv.bcf" | wc -l) invariant, $(bcftools view -H "${T}.snp.bcf" | wc -l) SNPs"
    rm -f "${T}".raw.bcf "${T}".inv.bcf* "${T}".snp.bcf*
    pixy --stats pi dxy fst --vcf "${T}.vcf.gz" --populations "${OUT}/populations.txt" \
        --window_size "$WIN" --chromosomes "$C" --interval_start "$S" --interval_end "$E" \
        --fst_type hudson --n_cores "$THREADS" --output_folder "${OUT}/pixy" --output_prefix "$TAG"
    if [[ "$C" == "chr_Z" ]]; then      # males only: no hemizygous females
        pixy --stats pi dxy fst --vcf "${T}.vcf.gz" --populations "${OUT}/populations_males.txt" \
            --window_size "$WIN" --chromosomes "$C" --interval_start "$S" --interval_end "$E" \
            --fst_type hudson --n_cores "$THREADS" --output_folder "${OUT}/pixy" --output_prefix "${TAG}_males"
    fi
    echo ">>> done $(date)"
    ;;
# -----------------------------------------------------------------------------
summary)
    need python3
    N=$(wc -l < "${OUT}/chunks.txt"); D=$(ls "${OUT}"/pixy/*_dxy.txt 2>/dev/null | grep -vc _males_ || true)
    echo ">>> pixy output for ${D} of ${N} chunks"
    python3 "${REPO}/07_pixy_summary.py" --pixy "${OUT}/pixy" --samples "$SAMPLES" --out "${OUT}/summary" \
        --min-sites "$MIN_SITES" --max-div "$MAX_DIV" --boot "$BOOT"
    # Figure 7: needs the reference-Z -> ptarmigan-Z alignment from step 06
    if cp "${PROJECT_DIR}/z_inversions/ref_to_ptarmigan.paf" "${OUT}/summary/" 2>/dev/null; then
        python3 "${REPO}/07_fig7_diversity_divergence.py" "${OUT}/summary"
    else
        echo "  (no z_inversions/ref_to_ptarmigan.paf: run step 06 with PTARMIGAN_FA for Figure 7)"
    fi
    ;;
*)
    echo "Usage: bash $0 samples | sbatch --export=ALL,STAGE=ref $0 | bash $0 map | bash $0 call | STAGE=summary sbatch --export=ALL $0" >&2
    exit 1 ;;
esac
