#!/bin/bash
# =============================================================================
# SLURM JOB SUBMISSION: PCA + RUNS OF HOMOZYGOSITY (combined cohort)
# Step 08 — requires 06_downsample_and_finalize.sh (final_cramlist.txt) and,
# for the heterozygosity aggregation step, 07_heterozygosity_array.sh to
# have completed for all samples. Input is CRAM (ANGSD's -bam flag accepts
# a list of CRAM paths the same way it does BAM — it reads them through
# the same htslib backend regardless of format).
#
# Replicates beagle.sh + pca.sh + ROH.sh + rohparser.py from
# https://github.com/Andrew-N-Black/LEPC-popgen, with two deliberate,
# explicitly-requested deviations from the original:
#   - Whole-genome analysis, not the original's 100kb-window reference
#     subset (no windowing scheme to replicate/fabricate).
#   - ANGSD -doGlf2 (beagle) is parallelized per chromosome instead of per
#     100kb window, since we're not subsetting the genome — this still
#     covers 100% of the genome, just chunked for tractability.
#
# Everything else (ANGSD/pcangsd/bcftools flags) matches the original
# scripts exactly, verified against their actual source rather than
# guessed — see the flag comments at each step below.
#
# -minInd is set to round(0.75 x N), matching the ~75% stringency the
# original used (-minInd 348 of their N=~464) rather than a hardcoded
# number, since N here depends on how many NEW samples end up sequenced.
#
# This version merges in the ROH pipeline fixes developed and confirmed
# working separately (originally split out as 09_roh.sh) and adds a third
# analysis, ROHan (Renaud et al. 2019), which estimates heterozygosity/ROH
# directly from BAM/CRAM via its own genotype-likelihood model rather than
# from called genotypes -- a useful cross-check against the ANGSD/bcftools
# approach above. ROHan runs last since its module environment (a clean
# gcc/samtools setup) is incompatible with the angsd/pcangsd/bcftools
# modules the earlier steps need, and there's no reason to juggle both at
# once when nothing later in the script needs the earlier modules again.
# install_rohan.sh must have been run on the login node first.
#
# USAGE:
#   sbatch 08_pca_roh.sh
# =============================================================================
#SBATCH --job-name=old.new_pca_roh
#SBATCH --output=logs/%x_%j.out
#SBATCH --error=logs/%x_%j.err
#SBATCH -A dewoody
#SBATCH -t 10-00:00:00
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=64
#SBATCH --mem=250G
#SBATCH -p cpu
#SBATCH --mail-type=BEGIN,END,FAIL
#SBATCH --mail-user=blackan@purdue.edu

set -euo pipefail

# =============================================================================
# USER-DEFINED VARIABLES
# =============================================================================
PROJECT_DIR="${CLUSTER_SCRATCH}/GROUSE/old_vs_new"
REF_FASTA="${PROJECT_DIR}/ref/GCF_026119805.1_pur_lepc_1.0_genomic.fna"
FINAL_CRAMLIST="${PROJECT_DIR}/final_cramlist.txt"
HET_DIR="${PROJECT_DIR}/heterozygosity"

BEAGLE_DIR="${PROJECT_DIR}/beagle"
PCA_DIR="${PROJECT_DIR}/pca"
ROH_DIR="${PROJECT_DIR}/roh"

# rohparser.py — vendored verbatim from the original repo rather than
# reimplemented, so ROH size-class/FROH logic matches exactly. Its one
# hardcoded path (a .fai file, for total genome length) is patched below
# to point at our reference instead of the original's (different cluster).
ROHPARSER_URL="https://raw.githubusercontent.com/Andrew-N-Black/LEPC-popgen/main/analysis/rohparser.py"
ROHPARSER="${ROH_DIR}/rohparser.py"
ROHPARSER_ORIG_FAI="${PROJECT_DIR}/ref/GCF_026119805.1_pur_lepc_1.0_genomic.fna.fai"

THREADS=$SLURM_CPUS_PER_TASK
ROH_PARALLEL_JOBS=8

# ROHan-specific config (Step 8, run last -- see module note there)
ROHAN_BIN="${PROJECT_DIR}/tools/ROHan/bin/rohan"
GSL_PREFIX="${PROJECT_DIR}/tools/gsl"
ROHAN_OUT_DIR="${PROJECT_DIR}/results/rohan"
ROHAN_THREADS=16
# ROHan's expected within-ROH heterozygosity rate parameter. ROHan's own
# examples use something on the order of 2e-5; adjust if your species'
# expected mutation rate differs substantially. Confirm the flag name is
# still --rohmu with `rohan --help` if ROHan is ever rebuilt/updated.
ROHMU=2e-5

mkdir -p logs "$BEAGLE_DIR" "$PCA_DIR" "$ROH_DIR" "$ROHAN_OUT_DIR"

echo ">>> 08_pca_roh.sh"
echo ">>> Start time: $(date)"

if [[ ! -f "$FINAL_CRAMLIST" ]]; then
    echo "ERROR: ${FINAL_CRAMLIST} not found. Run 06_downsample_and_finalize.sh first."
    exit 1
fi
if [[ ! -f "${REF_FASTA}.fai" ]]; then
    echo "ERROR: ${REF_FASTA}.fai not found. Run 05_combined_alignment_array.sh's prep step first."
    exit 1
fi

N_SAMPLES=$(wc -l < "$FINAL_CRAMLIST")
MININD=$(awk -v n="$N_SAMPLES" 'BEGIN { printf "%d", (n * 0.75) + 0.5 }')
echo ">>> N samples : ${N_SAMPLES}"
echo ">>> minInd    : ${MININD} (75% of N, matching the original's ~75% stringency)"

# =============================================================================
# ANGSD/pcangsd/bcftools ENVIRONMENT (Steps 0-6 below)
# =============================================================================
ml biocontainers
ml bcftools
ml angsd/0.940
ml pcangsd
ml htslib
# RCAC's xalt accounting hook injects LD_PRELOAD (libxalt_init.so) into
# every command, including containerized ones. singularity forwards it
# into the container by default, and the container's older glibc lacks
# the GLIBC_2.33/2.34 symbols that library needs, so angsd/pcangsd abort
# before running. Blanking it inside the container via these two env vars
# is the reliable fix (a plain `unset LD_PRELOAD` on the host doesn't
# hold — xalt is sticky and re-injects it):
export SINGULARITYENV_LD_PRELOAD=""
export APPTAINERENV_LD_PRELOAD=""

# =============================================================================
# STEP 0: Aggregate per-sample heterozygosity results from 07 (each array
# task there wrote its own file to avoid a shared-file race).
# =============================================================================
echo ">>> Step 0: Aggregating heterozygosity results"

HET_SUMMARY="${PROJECT_DIR}/heterozygosity_summary.tsv"
echo -e "sample_id\theterozygosity" > "$HET_SUMMARY"
find "$HET_DIR" -name "*_heterozygosity.txt" -exec cat {} + >> "$HET_SUMMARY" 2>/dev/null || true
echo "  Wrote ${HET_SUMMARY} ($(($(wc -l < "$HET_SUMMARY") - 1)) samples)"

# =============================================================================
# STEP 1: Genotype likelihoods (beagle format), parallelized per chromosome
# Flags match the original beagle.sh exactly (GL model, major/minor, MAF,
# quality, triallelic/SNP filtering) — only -minInd is recomputed for N,
# and region scope is per-chromosome instead of per-100kb-window.
# =============================================================================
echo ">>> Step 1: ANGSD genotype likelihoods (beagle format)"

CHROM_LIST="${BEAGLE_DIR}/chroms.txt"
cut -f1 "${REF_FASTA}.fai" > "$CHROM_LIST"
N_CHROMS=$(wc -l < "$CHROM_LIST")
echo "  ${N_CHROMS} chromosomes/contigs to process"

BEAGLE_THREADS_PER_JOB=$(( THREADS / ROH_PARALLEL_JOBS > 0 ? THREADS / ROH_PARALLEL_JOBS : 1 ))

run_beagle_chrom() {
    local chrom="$1"
    local out="${BEAGLE_DIR}/${chrom}"
    if [[ -f "${out}.beagle.gz" ]]; then
        return 0
    fi
    angsd -bam "$FINAL_CRAMLIST" -ref "$REF_FASTA" -r "${chrom}:" \
        -GL 1 -doGlf 2 -doMajorMinor 1 -doMaf 1 -minMaf 0.01 -minQ 30 \
        -skipTriallelic 1 -SNP_pval 1e-6 -minInd "$MININD" \
        -P "$BEAGLE_THREADS_PER_JOB" -out "$out"
}
# angsd is a bash function (from `ml angsd/0.940`'s Lmod setup, wrapping the
# singularity call), not a real binary on PATH. Functions don't propagate
# into the fresh bash process `xargs ... bash -c` spawns unless each one is
# individually exported with `export -f` — exporting run_beagle_chrom alone
# is not enough, since its body calls angsd, which was never exported.
export -f angsd
export -f run_beagle_chrom
export FINAL_CRAMLIST REF_FASTA BEAGLE_DIR MININD BEAGLE_THREADS_PER_JOB

xargs -a "$CHROM_LIST" -I{} -P "$ROH_PARALLEL_JOBS" bash -c 'run_beagle_chrom "$@"' _ {}

echo ">>> Step 1b: Concatenating per-chromosome beagle files"

FINAL_BEAGLE="${BEAGLE_DIR}/final.beagle.gz"
if [[ ! -f "$FINAL_BEAGLE" ]]; then
    FIRST=1
    > "${BEAGLE_DIR}/final.beagle"
    while IFS= read -r chrom; do
        f="${BEAGLE_DIR}/${chrom}.beagle.gz"
        [[ ! -f "$f" ]] && { echo "  WARNING: missing ${f} — skipping" >&2; continue; }
        if [[ "$FIRST" -eq 1 ]]; then
            zcat "$f" >> "${BEAGLE_DIR}/final.beagle"
            FIRST=0
        else
            zcat "$f" | tail -n +2 >> "${BEAGLE_DIR}/final.beagle"
        fi
    done < "$CHROM_LIST"
    gzip "${BEAGLE_DIR}/final.beagle"
fi

echo "  Final beagle file: ${FINAL_BEAGLE}"

# =============================================================================
# STEP 2: PCA + inbreeding (pcangsd) — flags match pca.sh exactly
# =============================================================================
echo ">>> Step 2: pcangsd"

pcangsd -b "$FINAL_BEAGLE" -o "${PCA_DIR}/final" --threads "$THREADS" --minMaf 0.01 --admix

pcangsd -b "$FINAL_BEAGLE" -o "${PCA_DIR}/final_inbreed" --threads "$THREADS" --minMaf 0.01 \
    --maf_tole 1e-9 --tole 1e-9 --inbreedSamples --inbreedSites \
    --iter 5000 --maf_iter 5000 --inbreed_iter 5000 --inbreed_tole 1e-9

echo "  PCA output      : ${PCA_DIR}/final.cov (+ .admix.Q etc.)"
echo "  Inbreeding output: ${PCA_DIR}/final_inbreed.*"

# =============================================================================
# STEP 3: ANGSD genome-wide variant calling -> BCF (flags match ROH.sh
# exactly, extracted directly from its source — no -doGeno needed)
# =============================================================================
echo ">>> Step 3: ANGSD variant calling (BCF output)"

JOINT_OUT="${ROH_DIR}/joint"
JOINT_BCF="${JOINT_OUT}.bcf"

if [[ ! -f "$JOINT_BCF" ]]; then
    angsd -bam "$FINAL_CRAMLIST" -ref "$REF_FASTA" \
        -GL 1 -dobcf 1 -dopost 1 -domajorminor 1 -domaf 1 \
        -minQ 30 -SNP_pval 1e-6 -P "$THREADS" -out "$JOINT_OUT"
else
    echo "  ${JOINT_BCF} already exists -- skipping ANGSD call. Delete it first for a clean rerun."
fi

if [[ ! -f "$JOINT_BCF" ]]; then
    echo "ERROR: ANGSD did not produce expected output: ${JOINT_BCF}"
    exit 1
fi

# =============================================================================
# STEP 4: Allele frequency file for bcftools roh
# =============================================================================
echo ">>> Step 4: Building allele-frequency file"

FREQS="${ROH_DIR}/freqs.tab.gz"
if [[ ! -f "$FREQS" ]]; then
    bcftools query -f '%CHROM\t%POS\t%REF,%ALT\t%AF\n' "$JOINT_BCF" | bgzip -c > "$FREQS"
    tabix -s1 -b2 -e2 "$FREQS"
else
    echo "  ${FREQS} already exists -- skipping."
fi

# =============================================================================
# STEP 5: bcftools roh (flags match ROH.sh exactly), streamed through grep
# so the multi-GB per-site ST output never touches disk -- only the RG
# (called-region) lines, which is all downstream parsing actually needs.
# =============================================================================
echo ">>> Step 5: bcftools roh"

ROH_RG_ONLY="${ROH_DIR}/ROH_GROUSE_PL_regions.txt"
bcftools roh --AF-file "$FREQS" --threads "$THREADS" "$JOINT_BCF" \
    | grep "^RG" > "$ROH_RG_ONLY"

echo "  RG (called-region) lines: ${ROH_RG_ONLY}"
echo "  $(wc -l < "$ROH_RG_ONLY") regions called across all samples"

# =============================================================================
# STEP 6: Per-sample ROH parsing with rohparser.py (vendored from the
# original repo, patched to use our reference's .fai for genome length)
# =============================================================================
echo ">>> Step 6: Per-sample ROH parsing"

if [[ ! -f "$ROHPARSER" ]]; then
    echo ">>> Downloading rohparser.py"
    wget -q -O "$ROHPARSER" "$ROHPARSER_URL"
    sed -i "s|${ROHPARSER_ORIG_FAI}|${REF_FASTA}.fai|g" "$ROHPARSER"
    # sed doesn't error or warn if ROHPARSER_ORIG_FAI didn't actually match
    # anything in the downloaded file -- it just silently leaves the
    # original (wrong-cluster) path in place. Fail loudly instead of
    # discovering this later as silently-wrong ROH results.
    if ! grep -qF "${REF_FASTA}.fai" "$ROHPARSER"; then
        echo "ERROR: rohparser.py patch did not take -- ROHPARSER_ORIG_FAI" >&2
        echo "  ('${ROHPARSER_ORIG_FAI}') was not found verbatim in the" >&2
        echo "  downloaded script. Check the source hasn't changed its" >&2
        echo "  hardcoded path, update ROHPARSER_ORIG_FAI to match, and" >&2
        echo "  delete ${ROHPARSER} to force a fresh download+patch." >&2
        exit 1
    fi
fi

# Single pass over the RG-only file, splitting by sample. NOTE: bcftools
# roh's sample column here is whatever the BCF's own header used as the
# sample name -- and ANGSD's -dobcf output uses the full BAM/CRAM file path
# as that identifier, not a bare sample ID (confirmed from the actual RG
# lines: column 2 is a full "/scratch/.../crams/F10.md.dedup_q20.cram"
# path). Using that path verbatim as a filename produces a broken, doubled
# path, so derive a clean sample ID from its basename instead, stripping
# this project's established CRAM suffix.
awk -v dir="$ROH_DIR" '
{
    n = split($2, parts, "/")
    sample = parts[n]
    gsub(/\.md\.dedup_q20\.cram$/, "", sample)
    print > (dir"/"sample"ROH.txt")
}' "$ROH_RG_ONLY"

N_SAMPLE_FILES=$(find "$ROH_DIR" -maxdepth 1 -name "*ROH.txt" ! -empty | wc -l)
echo "  Split into ${N_SAMPLE_FILES} non-empty per-sample files (expected ${N_SAMPLES})"
if [[ "$N_SAMPLE_FILES" -ne "$N_SAMPLES" ]]; then
    echo "  WARNING: sample-file count doesn't match N_SAMPLES -- a sample may" >&2
    echo "  have zero called ROH regions (possible, not necessarily a bug)," >&2
    echo "  or something upstream is off. Compare against:" >&2
    echo "    bcftools query -l ${JOINT_BCF}" >&2
fi

run_rohparser() {
    # rohparser.py builds its own input path internally from a hardcoded
    # directory + a bare filename (matching its documented usage: `cd` into
    # the ROH directory, then `python ROHparser.py SAMPLEROH.txt`). Passing
    # it a full path instead -- as a naive `find`-based invocation would --
    # makes it concatenate a doubled, nonexistent path and fail. So: cd into
    # the file's directory and pass only the basename.
    local roh_file="$1"
    local bn
    bn=$(basename "$roh_file")
    (cd "$(dirname "$roh_file")" && python3 "$ROHPARSER" "$bn") > "${roh_file}_results.txt"
}
export -f run_rohparser
export ROHPARSER

find "$ROH_DIR" -maxdepth 1 -name "*ROH.txt" ! -empty \
    | xargs -I{} -P "$ROH_PARALLEL_JOBS" bash -c 'run_rohparser "$@"' _ {}

N_ROH_RESULTS=$(find "$ROH_DIR" -maxdepth 1 -name "*ROH.txt_results.txt" | wc -l)
echo "  Parsed ANGSD/bcftools ROH results for ${N_ROH_RESULTS} samples"

# =============================================================================
# STEP 7: ROHan (Renaud et al. 2019) -- per-sample heterozygosity/ROH
# estimated directly from BAM/CRAM via its own genotype-likelihood model,
# as a cross-check against the ANGSD/bcftools estimates above. Does NOT
# take a VCF/BCF as input by design (that's the point of the tool -- it
# stays upstream of hard-called genotypes). Requires install_rohan.sh to
# have been run on the login node first.
#
# Runs LAST and does its own `module --force purge`: its build needs a
# plain gcc/samtools environment incompatible with the angsd/pcangsd/
# bcftools modules loaded above, and nothing after this point needs those
# modules again.
# =============================================================================
echo ">>> Step 7: ROHan per-sample analysis"

if [ ! -x "$ROHAN_BIN" ]; then
    echo "ERROR: ROHan binary not found at $ROHAN_BIN" >&2
    echo "  Run install_rohan.sh on the login node first." >&2
    exit 1
fi

module --force purge
module load gcc/14.1.0
module load biocontainers
module load samtools

# rohan was linked against a custom-built GSL (no 'gsl' module exists on
# Gautschi), so it needs this to find libgsl.so at runtime, not just at
# build time -- see install_rohan.sh.
export LD_LIBRARY_PATH="$GSL_PREFIX/lib${LD_LIBRARY_PATH:+:$LD_LIBRARY_PATH}"

while IFS= read -r CRAM; do
    SAMPLE=$(basename "$CRAM" | sed -E 's/\.md\.dedup_q20\.cram$//')
    OUT_PREFIX="${ROHAN_OUT_DIR}/${SAMPLE}"
    if [[ -f "${OUT_PREFIX}.hEst" ]]; then
        echo "  ${SAMPLE}: ${OUT_PREFIX}.hEst already exists -- skipping."
        continue
    fi
    echo "=== ROHan: $SAMPLE ==="

    # ROHan's most reliably supported input format is BAM. Rather than rely
    # on ROHan's own CRAM/reference handling (undocumented in what we could
    # verify), convert to a temporary indexed BAM first -- slower, but
    # removes any ambiguity about reference resolution.
    TMP_BAM="${ROHAN_OUT_DIR}/${SAMPLE}.tmp.bam"
    samtools view -@ "$ROHAN_THREADS" -b -T "$REF_FASTA" -o "$TMP_BAM" "$CRAM"
    samtools index "$TMP_BAM"

    "$ROHAN_BIN" \
        -t "$ROHAN_THREADS" \
        --rohmu "$ROHMU" \
        -o "$OUT_PREFIX" \
        "$REF_FASTA" "$TMP_BAM"

    rm -f "$TMP_BAM" "${TMP_BAM}.bai"
    echo "  Done: ${OUT_PREFIX}.*"
done < "$FINAL_CRAMLIST"

echo ""
echo ">>> PCA + ROH + ROHan analysis complete."
echo "    Heterozygosity      : ${HET_SUMMARY}"
echo "    PCA                 : ${PCA_DIR}/final.cov"
echo "    Inbreeding          : ${PCA_DIR}/final_inbreed.*"
echo "    Joint BCF           : ${JOINT_BCF}"
echo "    ANGSD/bcftools ROH  : ${ROH_DIR}/*ROH.txt_results.txt"
echo "    ROHan (per sample)  : ${ROHAN_OUT_DIR}/<sample>.hEst, .mid.hmmp, .mid.ROH"
echo ">>> End time: $(date)"
