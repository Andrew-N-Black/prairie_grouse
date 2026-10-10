#!/bin/bash
# =============================================================================
# SLURM ARRAY JOB: DE NOVO ASSEMBLY AND Hi-C SCAFFOLDING — hifiasm + yahs
# Step 02 — the reference-free half of the pipeline. One array task per sample;
# n=23 samples across 3 species (STGR, LEPC, GRPC).
#
# This script stops at the yahs scaffolds deliberately. Everything downstream of
# here is reference-guided, and the reference is a choice that can be revisited
# without re-running any of this: step 05 orders these same scaffolds into
# pseudo-chromosomes. Keeping the split here means a change of reference costs a
# day, not a month.
#
# Per sample, per haplotype (hap1/hap2):
#   1.  HiFiAdapterFilt   — adapter filtering of raw HiFi BAMs (per sample)
#   2.  hifiasm           — Hi-C-phased assembly, +ONT UL when available
#   3.  GFA -> FASTA
#   3b. QUAST             — post-hifiasm contig stats
#   4.  bwa mem -5SP      — Hi-C reads -> contigs (MAPQ >= 20)
#   5.  yahs              — Hi-C scaffolding
#   5b. QUAST             — post-yahs scaffold stats
#   6.  juicer pre + juicer_tools — .hic/.assembly for Juicebox curation
#
# OUTPUT consumed by step 05:
#   ${PROJECT_DIR}/yahs/<SAMPLE>.<HAP>_scaffolds_final.fa
#
# Every step is resume-safe: completed outputs are detected and skipped, so a
# failed task can simply be resubmitted.
#
# INPUT MANIFEST (tab-separated, header required, see assembly_manifest.tsv):
#   sample_id  species  hifi_bams  hic_r1  hic_r2  ont_ul
#   - hifi_bams : raw, unaligned PacBio HiFi BAM(s) (*.hifi_reads.bc####.bam),
#                 comma-separated if multiple SMRT cells
#   - hic_r1/r2 : paired-end Hi-C fastq.gz
#   - ont_ul    : ONT ultra-long fastq.gz, or "NA"
#
# USAGE:
#   N=$(grep -v '^#' assembly_manifest.tsv | tail -n +2 | grep -c .)
#   sbatch --array=0-$((N-1))%6 02_assembly_scaffold_array.sh
#   (%6 caps concurrent tasks — tune to fair-share/node availability)
#
#   Rerun a single failed task, e.g. index 7:
#   sbatch --array=7 02_assembly_scaffold_array.sh
# =============================================================================
#SBATCH --job-name=grouse_asm
#SBATCH --output=logs/%x_%A_%a.out
#SBATCH --error=logs/%x_%A_%a.err
#SBATCH -A fnrdewoody
#SBATCH -t 10-00:00:00
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=48
# 160G: measured peak RSS across 17 array tasks was 128-143 GB (sacct MaxRSS,
# whole-cgroup, at 128 threads with unpinned samtools sort). With 48 threads
# and `samtools sort -m 1G` below, the sort contribution is capped near 48 GB,
# so the real peak should land comfortably under this. On a 377 GB node this
# fits 2 tasks concurrently. Re-check `sacct -o MaxRSS` after a few samples;
# if the peak drops below ~120 GB, 3 tasks per node becomes viable.
#SBATCH --mem=160G
#SBATCH -p cpu
#SBATCH --mail-type=BEGIN,END,FAIL
#SBATCH --mail-user=blackan@purdue.edu

# =============================================================================
# ENVIRONMENT SETUP
# =============================================================================
set -euo pipefail

# Defensive: if the submitting shell had anaconda loaded, --export=ALL (the
# sbatch default) carries that Lmod state into the job, and liftoff will
# refuse to load later. Start from a known state.
module unload anaconda 2>/dev/null || true

ml biocontainers
ml quast
ml hifiasm
ml bwa
ml samtools/1.22.1
unset LD_PRELOAD
# unset LD_PRELOAD: RCAC's XALT usage-tracking library is injected via
# LD_PRELOAD and fails on some nodes (GLIBC_2.33/2.34 mismatch), which can
# kill subshells under `set -e`. It's accounting only — safe to drop.
#
# =============================================================================
# USER-DEFINED VARIABLES
# =============================================================================
MANIFEST="${SLURM_SUBMIT_DIR}/assembly_manifest.tsv"

# Shared across every step of the pipeline. The conda environments built here
# are reused by steps 05 and 09, so this is the only place they are created.
PROJECT_DIR="${CLUSTER_SCRATCH}/GROUSE/grouse_asm"
REF_DIR="${PROJECT_DIR}/ref"

FILT_DIR="${PROJECT_DIR}/hifi_filtered"
ASM_DIR="${PROJECT_DIR}/hifiasm"
SCAFFOLD_DIR="${PROJECT_DIR}/yahs"
QC_DIR="${PROJECT_DIR}/qc"

HIC_MIN_MAPQ=20

CONDA_ENVS_DIR="${PROJECT_DIR}/conda_envs"

YAHS_VERSION="1.2.2"
YAHS_ENV_DIR="${CONDA_ENVS_DIR}/yahs-${YAHS_VERSION}"
YAHS_BIN="${YAHS_ENV_DIR}/bin/yahs"
JUICER_BIN="${YAHS_ENV_DIR}/bin/juicer"   # yahs's own JBAT pre-processor

# juicer_tools via bioconda, pinned. Manual jar downloads proved unreliable
# (the classic S3 mirror 403s) and v3.0.0 dropped the classic `pre` syntax.
JUICER_TOOLS_VERSION="2.20.00"
JUICER_TOOLS_ENV_DIR="${CONDA_ENVS_DIR}/juicertools-${JUICER_TOOLS_VERSION}"

PRETEXT_ENV_DIR="${CONDA_ENVS_DIR}/pretext"
PRETEXTMAP_BIN="${PRETEXT_ENV_DIR}/bin/PretextMap"
PRETEXTSNAPSHOT_BIN="${PRETEXT_ENV_DIR}/bin/PretextSnapshot"

DOTPLOT_ENV_DIR="${CONDA_ENVS_DIR}/dotplot-python"
DOTPLOT_PYTHON_BIN="${DOTPLOT_ENV_DIR}/bin/python3"

# ACTIVATED (not absolute-path invoked) in Step 1 — its bioconda packaging
# relies on conda's activate hook to put its adapter database on PATH.
HIFIADAPTERFILT_ENV_DIR="${CONDA_ENVS_DIR}/hifiadapterfilt"
HIFIADAPTERFILT_SCRIPT="${HIFIADAPTERFILT_ENV_DIR}/bin/hifiadapterfilt.sh"

TIDK_ENV_DIR="${CONDA_ENVS_DIR}/tidk"
TIDK_BIN="${TIDK_ENV_DIR}/bin/tidk"

THREADS=$SLURM_CPUS_PER_TASK

mkdir -p logs "$CONDA_ENVS_DIR" "$REF_DIR" "$FILT_DIR" "$ASM_DIR" "$SCAFFOLD_DIR" \
    "$QC_DIR" "${QC_DIR}/quast"

# =============================================================================
# SHARED ONE-TIME SETUP (conda envs)
# Serialized with flock: with many array tasks starting at once, two tasks
# running `conda create` into the same prefix (or writing the same map file)
# simultaneously would corrupt it. The first task to get the lock does the
# work; the rest wait, then find everything already present and skip.
# =============================================================================
exec 9>"${CONDA_ENVS_DIR}/.setup.lock"
flock 9

NEED_ANACONDA=false
[[ ! -x "$YAHS_BIN" || ! -x "$JUICER_BIN" ]]              && NEED_ANACONDA=true
[[ ! -x "$PRETEXTMAP_BIN" || ! -x "$PRETEXTSNAPSHOT_BIN" ]] && NEED_ANACONDA=true
[[ ! -x "$DOTPLOT_PYTHON_BIN" ]]                          && NEED_ANACONDA=true
[[ ! -f "$HIFIADAPTERFILT_SCRIPT" ]]                      && NEED_ANACONDA=true
[[ ! -d "$JUICER_TOOLS_ENV_DIR" ]]                        && NEED_ANACONDA=true
[[ ! -x "$TIDK_BIN" ]]                                    && NEED_ANACONDA=true

if [[ "$NEED_ANACONDA" == true ]]; then
    ml anaconda/2025.12-py313

    if [[ ! -x "$YAHS_BIN" || ! -x "$JUICER_BIN" ]]; then
        echo ">>> Installing yahs v${YAHS_VERSION}"
        conda create --yes --override-channels --prefix "$YAHS_ENV_DIR" -c bioconda -c conda-forge "yahs=${YAHS_VERSION}"
    fi
    if [[ ! -x "$PRETEXTMAP_BIN" || ! -x "$PRETEXTSNAPSHOT_BIN" ]]; then
        echo ">>> Installing PretextMap/PretextSnapshot"
        conda create --yes --override-channels --prefix "$PRETEXT_ENV_DIR" -c bioconda -c conda-forge \
            pretextmap=0.2.4 pretextsnapshot=0.0.7
    fi
    if [[ ! -x "$DOTPLOT_PYTHON_BIN" ]]; then
        echo ">>> Installing Python 3 + matplotlib"
        conda create --yes --override-channels --prefix "$DOTPLOT_ENV_DIR" -c conda-forge python=3.11 matplotlib
    fi
    if [[ ! -f "$HIFIADAPTERFILT_SCRIPT" ]]; then
        echo ">>> Installing HiFiAdapterFilt"
        conda create --yes --override-channels --prefix "$HIFIADAPTERFILT_ENV_DIR" -c bioconda -c conda-forge hifiadapterfilt
    fi
    if [[ ! -d "$JUICER_TOOLS_ENV_DIR" ]]; then
        echo ">>> Installing juicertools v${JUICER_TOOLS_VERSION}"
        conda create --yes --override-channels --prefix "$JUICER_TOOLS_ENV_DIR" -c bioconda -c conda-forge "juicertools=${JUICER_TOOLS_VERSION}"
    fi
    if [[ ! -x "$TIDK_BIN" ]]; then
        echo ">>> Installing tidk"
        conda create --yes --override-channels --prefix "$TIDK_ENV_DIR" -c bioconda -c conda-forge tidk
    fi

    module unload anaconda/2025.12-py313
fi

flock -u 9
exec 9>&-

# Verify every tool is actually in place.
for BIN in "$YAHS_BIN" "$JUICER_BIN" "$PRETEXTMAP_BIN" "$PRETEXTSNAPSHOT_BIN" "$DOTPLOT_PYTHON_BIN" "$TIDK_BIN"; do
    if [[ ! -x "$BIN" ]]; then
        echo "ERROR: expected tool not found/executable: ${BIN}"
        exit 1
    fi
done
if [[ ! -f "$HIFIADAPTERFILT_SCRIPT" ]]; then
    echo "ERROR: HiFiAdapterFilt not found: ${HIFIADAPTERFILT_SCRIPT}"
    exit 1
fi

# Resolve juicer_tools: bioconda wrapper script if present, else the jar.
if [[ -x "${JUICER_TOOLS_ENV_DIR}/bin/juicer_tools" ]]; then
    JUICER_TOOLS_CMD=("${JUICER_TOOLS_ENV_DIR}/bin/juicer_tools" -Xmx32G)
else
    JUICER_TOOLS_JAR_FOUND=$(find "$JUICER_TOOLS_ENV_DIR" -iname "juicer_tools*.jar" 2>/dev/null | head -n1)
    if [[ -n "$JUICER_TOOLS_JAR_FOUND" ]]; then
        JUICER_TOOLS_CMD=(java -Xmx32G -jar "$JUICER_TOOLS_JAR_FOUND")
    else
        echo "ERROR: no juicer_tools executable or jar under ${JUICER_TOOLS_ENV_DIR}"
        exit 1
    fi
fi

# =============================================================================
# RESOLVE SAMPLE FOR THIS ARRAY TASK
# =============================================================================
if [[ -z "${SLURM_ARRAY_TASK_ID:-}" ]]; then
    echo "ERROR: SLURM_ARRAY_TASK_ID is not set."
    echo "Submit with: sbatch --array=0-N 02_assembly_scaffold_array.sh"
    exit 1
fi
if [[ ! -f "$MANIFEST" ]]; then
    echo "ERROR: Manifest not found: ${MANIFEST}"
    exit 1
fi

mapfile -t ROWS < <(grep -v '^#' "$MANIFEST" | tail -n +2 | grep -v '^[[:space:]]*$')
LINE="${ROWS[$SLURM_ARRAY_TASK_ID]:-}"
if [[ -z "$LINE" ]]; then
    echo "ERROR: No manifest row at index ${SLURM_ARRAY_TASK_ID} (${#ROWS[@]} samples in manifest)"
    exit 1
fi

IFS=$'\t' read -r SAMPLE SPECIES HIFI_BAMS HIC_R1 HIC_R2 ONT_UL <<< "$LINE"
ONT_UL="${ONT_UL:-NA}"
ONT_UL="${ONT_UL%$'\r'}"   # tolerate Windows line endings in the manifest

if [[ -z "$SAMPLE" || -z "$SPECIES" || -z "$HIFI_BAMS" || -z "$HIC_R1" || -z "$HIC_R2" ]]; then
    echo "ERROR: Malformed manifest row at index ${SLURM_ARRAY_TASK_ID}: ${LINE}"
    exit 1
fi

echo ">>> Array task ${SLURM_ARRAY_TASK_ID} -> sample: ${SAMPLE} (${SPECIES})"
echo ">>> HiFi BAMs  : ${HIFI_BAMS}"
echo ">>> Hi-C reads : ${HIC_R1} / ${HIC_R2}"
echo ">>> ONT UL     : ${ONT_UL}"
echo ">>> Running on : $(hostname)"
echo ">>> CPUs       : ${THREADS}"
echo ">>> juicer_tools: ${JUICER_TOOLS_CMD[*]}"
echo ">>> Start time : $(date)"

IFS=',' read -ra HIFI_BAM_ARR <<< "$HIFI_BAMS"
for F in "${HIFI_BAM_ARR[@]}" "$HIC_R1" "$HIC_R2"; do
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
# STEP 1: HiFiAdapterFilt (per sample)
#
# IMPORTANT (learned the hard way): hifiadapterfilt.sh sets `outdir=$(pwd)` by
# default and uses TWO different variables internally — it converts the BAM to
# FASTQ/FASTA next to the BAM (its `read_path_str`), but hands BLAST a path
# built from `outdir`. If those two directories differ, BLAST silently finds
# no query file, no blocklist is produced, and the "filtered" output is a
# VERBATIM COPY of the input with nothing removed. That failure is silent:
# the .filt.fastq.gz exists, so any existence-based guard happily skips it.
#
# Fix: symlink the BAM into the output directory and run entirely inside it
# (cd there, -o "."), so read_path_str and outdir are the same directory.
#
# NOTE ON DISK: this step converts each BAM to uncompressed FASTQ *and* FASTA.
# With 130-190 GB BAMs that is roughly 400-600 GB of transient scratch per
# sample, on top of the BAM itself. The script removes the FASTA itself; the
# FASTQ and the symlink are cleaned up below once filtering succeeds.
# =============================================================================
echo ">>> Step 1: Adapter filtering raw HiFi BAMs (HiFiAdapterFilt)"

FILT_SAMPLE_DIR="${FILT_DIR}/${SAMPLE}"
mkdir -p "$FILT_SAMPLE_DIR"

FILT_HIFI_ARR=()
NEED_FILTER=false
for RAW_HIFI_BAM in "${HIFI_BAM_ARR[@]}"; do
    PREFIX=$(basename "$RAW_HIFI_BAM" .bam)
    [[ -f "${FILT_SAMPLE_DIR}/${PREFIX}.filt.fastq.gz" ]] || NEED_FILTER=true
done

if [[ "$NEED_FILTER" == true ]]; then
    ml anaconda/2025.12-py313
    conda activate "$HIFIADAPTERFILT_ENV_DIR"

    # --- BLAST database path shim (bioconda packaging bug) ---
    # hifiadapterfilt.sh line 7 derives its BLAST database location with:
    #   DBpath=$(echo $PATH | sed 's/:/\n/g' | grep "HiFiAdapterFilt/DB" | head -n 1)
    # i.e. it greps PATH for the literal string "HiFiAdapterFilt/DB", which is
    # the upstream GitHub layout. bioconda installs the database at
    # <env>/bin/DB (lowercase package name, "bin" not the package name), so
    # that grep matches NOTHING, DBpath ends up empty, and every blastn call
    # becomes `-db /pacbio_vectors_db` -> "BLAST Database error: No alias or
    # index file found". BLAST writes nothing, the blocklist is empty, and
    # NOTHING is filtered -- silently, while the stats file cheerfully reports
    # "0 adapter contaminated ccs reads (0% of total)".
    #
    # DBpath is assigned with a plain `=`, so exporting it here would just be
    # overwritten. Instead give that grep something to match: a symlink whose
    # path literally contains "HiFiAdapterFilt/DB", appended to PATH.
    HIFIADAPTERFILT_DB_SHIM="${CONDA_ENVS_DIR}/HiFiAdapterFilt"
    mkdir -p "$HIFIADAPTERFILT_DB_SHIM"
    ln -sfn "${HIFIADAPTERFILT_ENV_DIR}/bin/DB" "${HIFIADAPTERFILT_DB_SHIM}/DB"
    export PATH="${PATH}:${HIFIADAPTERFILT_DB_SHIM}/DB"

    if ! echo "$PATH" | tr ':' '\n' | grep -q "HiFiAdapterFilt/DB"; then
        echo "ERROR: BLAST database shim not on PATH — HiFiAdapterFilt would"
        echo "silently filter nothing. Expected: ${HIFIADAPTERFILT_DB_SHIM}/DB"
        exit 1
    fi
    if [[ ! -f "${HIFIADAPTERFILT_DB_SHIM}/DB/pacbio_vectors_db.nin" ]]; then
        echo "ERROR: BLAST database not found via shim:"
        echo "  ${HIFIADAPTERFILT_DB_SHIM}/DB/pacbio_vectors_db.nin"
        exit 1
    fi
fi

for RAW_HIFI_BAM in "${HIFI_BAM_ARR[@]}"; do
    PREFIX=$(basename "$RAW_HIFI_BAM" .bam)
    FILT_FASTQ="${FILT_SAMPLE_DIR}/${PREFIX}.filt.fastq.gz"
    BLOCKLIST="${FILT_SAMPLE_DIR}/${PREFIX}.blocklist"
    STATS="${FILT_SAMPLE_DIR}/${PREFIX}.stats"

    if [[ -f "$FILT_FASTQ" ]]; then
        echo "  ${PREFIX}: already filtered — skipping"
    else
        echo "  Filtering ${PREFIX}"
        ln -sf "$RAW_HIFI_BAM" "${FILT_SAMPLE_DIR}/${PREFIX}.bam"
        # Tee stderr to a file so BLAST failures can be detected: the script
        # backgrounds its blastn calls and never checks their exit status, so
        # a database error otherwise just yields an empty blocklist and a
        # clean-looking "0% contaminated" result.
        FILT_STDERR="${FILT_SAMPLE_DIR}/${PREFIX}.hifiadapterfilt.stderr"
        # BLAST_THREADS is deliberately NOT $THREADS. hifiadapterfilt passes -t
        # straight to `blastn -num_threads`, and at high thread counts BLAST
        # hits "CThread::Run() -- error creating thread" and dies PARTWAY
        # through the search. That leaves a truncated .contaminant.blastout, an
        # empty blocklist, and a stats file reporting a confident but wrong
        # "0% contaminated". BLAST also scales poorly past ~8-16 threads, so
        # capping it costs little. The BAM->FASTQ conversion is single-threaded
        # regardless and dominates this step's runtime.
        BLAST_THREADS=8
        (cd "$FILT_SAMPLE_DIR" && hifiadapterfilt.sh -p "$PREFIX" -o "." -t "$BLAST_THREADS") \
            2> >(tee "$FILT_STDERR" >&2)

        # Any NCBI/BLAST failure means the adapter search did not complete, so
        # the "contaminated reads" count cannot be trusted. Fail loudly rather
        # than assembling from reads that were never fully screened.
        if grep -qi "BLAST Database error\|No alias or index file found\|error creating thread\|NCBI C++ Exception" "$FILT_STDERR"; then
            echo "ERROR: BLAST did not complete — adapter screening is incomplete."
            echo "The read counts in ${PREFIX}.stats are NOT trustworthy."
            echo "See ${FILT_STDERR}"
            exit 1
        fi
    fi

    if [[ ! -f "$FILT_FASTQ" ]]; then
        echo "ERROR: HiFiAdapterFilt did not produce expected output: ${FILT_FASTQ}"
        exit 1
    fi

    # Hard validation: the silent-no-op failure mode above produces a
    # .filt.fastq.gz with NO blocklist and NO stats file. Refuse to continue
    # on an unfiltered copy rather than assembling from it.
    if [[ ! -f "$BLOCKLIST" || ! -s "$STATS" ]]; then
        echo "ERROR: HiFiAdapterFilt produced ${FILT_FASTQ} but no blocklist/stats."
        echo "That means BLAST never ran and NOTHING was filtered — the output is"
        echo "an unfiltered copy of the input. Not continuing."
        echo "Expected: ${BLOCKLIST} and ${STATS}"
        exit 1
    fi
    echo "  --- ${PREFIX} filtering summary ---"
    grep -E "Number of (ccs reads|adapter contaminated|ccs reads retained)" "$STATS" || cat "$STATS"

    # Drop the huge uncompressed intermediates and the BAM symlink.
    rm -f "${FILT_SAMPLE_DIR}/${PREFIX}.fastq" "${FILT_SAMPLE_DIR}/${PREFIX}.fasta" \
          "${FILT_SAMPLE_DIR}/${PREFIX}.fq" "${FILT_SAMPLE_DIR}/${PREFIX}.bam"

    FILT_HIFI_ARR+=("$FILT_FASTQ")
done

if [[ "$NEED_FILTER" == true ]]; then
    conda deactivate
    module unload anaconda/2025.12-py313
fi

echo "  Filtered HiFi reads: ${FILT_HIFI_ARR[*]}"

# =============================================================================
# STEP 2: hifiasm (per sample)
# =============================================================================
echo ">>> Step 2: hifiasm assembly"

OUT_PREFIX="${ASM_DIR}/${SAMPLE}"

if [[ -f "${OUT_PREFIX}.hic.hap1.p_ctg.gfa" && -f "${OUT_PREFIX}.hic.hap2.p_ctg.gfa" ]]; then
    echo "  hifiasm hap1/hap2 GFAs already exist — skipping"
else
    HIFIASM_CMD=(hifiasm -o "$OUT_PREFIX" -t "$THREADS" --h1 "$HIC_R1" --h2 "$HIC_R2")
    if [[ "$HAS_UL" == true ]]; then
        echo "  Ultra-long ONT reads detected — adding --ul"
        HIFIASM_CMD+=(--ul "$ONT_UL")
    fi
    HIFIASM_CMD+=("${FILT_HIFI_ARR[@]}")
    echo "  ${HIFIASM_CMD[*]}"
    "${HIFIASM_CMD[@]}"
fi

# =============================================================================
# STEPS 3-11: per haplotype
# =============================================================================
for HAP in hap1 hap2; do
    echo ""
    echo ">>> ===== ${SAMPLE} ${HAP} ====="

    # ---------------------------------------------------------------- Step 3
    echo ">>> Step 3 (${HAP}): GFA -> FASTA"
    GFA="${OUT_PREFIX}.hic.${HAP}.p_ctg.gfa"
    CONTIGS="${ASM_DIR}/${SAMPLE}.${HAP}.contigs.fa"

    if [[ ! -f "$GFA" ]]; then
        echo "ERROR: Expected hifiasm output not found: ${GFA}"
        exit 1
    fi
    if [[ -f "$CONTIGS" && -f "${CONTIGS}.fai" ]]; then
        echo "  contigs FASTA already exists — skipping"
    else
        awk '/^S/{print ">"$2"\n"$3}' "$GFA" > "$CONTIGS"
        samtools faidx "$CONTIGS"
    fi

    # --------------------------------------------------------------- Step 3b
    echo ">>> Step 3b (${HAP}): QUAST — post-hifiasm contigs"
    QUAST_HIFIASM_OUT="${QC_DIR}/quast/${SAMPLE}.${HAP}.post_hifiasm"
    if [[ -f "${QUAST_HIFIASM_OUT}/report.txt" ]]; then
        echo "  already exists — skipping"
    else
        quast.py -o "$QUAST_HIFIASM_OUT" -t "$THREADS" --large "$CONTIGS"
    fi

    # ---------------------------------------------------------------- Step 4
    echo ">>> Step 4 (${HAP}): Align Hi-C reads to contigs (MAPQ >= ${HIC_MIN_MAPQ})"
    if [[ -f "${CONTIGS}.bwt" ]]; then
        echo "  bwa index already exists — skipping"
    else
        bwa index "$CONTIGS"
    fi

    # Stays BAM: yahs reads it directly and has no CRAM support.
    HIC_BAM="${SCAFFOLD_DIR}/${SAMPLE}.${HAP}.hic2contigs.bam"
    if [[ -f "$HIC_BAM" ]]; then
        echo "  Hi-C BAM already exists — skipping"
    else
        bwa mem -5SP -t "$THREADS" "$CONTIGS" "$HIC_R1" "$HIC_R2" \
            | samtools view -@ "$THREADS" -buS -q "$HIC_MIN_MAPQ" - \
            | samtools sort -@ "$THREADS" -m 1G -n -o "${HIC_BAM}.part" -
        mv "${HIC_BAM}.part" "$HIC_BAM"
    fi

    # ---------------------------------------------------------------- Step 5
    echo ">>> Step 5 (${HAP}): yahs Hi-C scaffolding"
    YAHS_PREFIX="${SCAFFOLD_DIR}/${SAMPLE}.${HAP}"
    SCAFFOLDS="${YAHS_PREFIX}_scaffolds_final.fa"
    if [[ -f "$SCAFFOLDS" ]]; then
        echo "  yahs scaffolds already exist — skipping"
    else
        "$YAHS_BIN" "$CONTIGS" "$HIC_BAM" -o "$YAHS_PREFIX"
    fi
    if [[ ! -f "$SCAFFOLDS" ]]; then
        echo "ERROR: yahs did not produce expected output: ${SCAFFOLDS}"
        exit 1
    fi

    # --------------------------------------------------------------- Step 5b
    echo ">>> Step 5b (${HAP}): QUAST — post-yahs scaffolds"
    QUAST_YAHS_OUT="${QC_DIR}/quast/${SAMPLE}.${HAP}.post_yahs"
    if [[ -f "${QUAST_YAHS_OUT}/report.txt" ]]; then
        echo "  already exists — skipping"
    else
        quast.py -o "$QUAST_YAHS_OUT" -t "$THREADS" --large "$SCAFFOLDS"
    fi

    # ---------------------------------------------------------------- Step 6
    echo ">>> Step 6 (${HAP}): Juicebox .hic/.assembly from raw yahs scaffolding"
    JBAT_PREFIX="${SCAFFOLD_DIR}/${SAMPLE}.${HAP}_JBAT"
    if [[ -f "${JBAT_PREFIX}.hic" && -f "${JBAT_PREFIX}.assembly" ]]; then
        echo "  Juicebox files already exist — skipping"
    else
        "$JUICER_BIN" pre -a -o "$JBAT_PREFIX" \
            "${YAHS_PREFIX}.bin" "${YAHS_PREFIX}_scaffolds_final.agp" "${CONTIGS}.fai" \
            > "${JBAT_PREFIX}.log" 2>&1

        # -n skips normalization-vector calculation. juicer_tools 2.x throws a
        # NullPointerException in that stage (AddNorm/NormalizationCalculations)
        # AFTER the .hic body is already written, and JBAT curation does not
        # need those vectors. Failure here is non-fatal: this file is only for
        # optional manual curation, so it must not abort the whole sample.
        "${JUICER_TOOLS_CMD[@]}" pre -n \
            "${JBAT_PREFIX}.txt" "${JBAT_PREFIX}.hic.part" \
            <(grep PRE_C_SIZE "${JBAT_PREFIX}.log" | awk '{print $2" "$3}') \
            || echo "  WARNING: juicer_tools pre failed — continuing without the Juicebox .hic"

        if [[ -s "${JBAT_PREFIX}.hic.part" ]]; then
            mv "${JBAT_PREFIX}.hic.part" "${JBAT_PREFIX}.hic"
        else
            echo "  WARNING: no .hic produced for ${SAMPLE} ${HAP}; JBAT curation unavailable"
        fi
    fi
    echo "  ${JBAT_PREFIX}.hic + ${JBAT_PREFIX}.assembly (open together in Juicebox)"
done

echo ""
echo "============================================================"
echo ">>> Sample ${SAMPLE} (${SPECIES}) complete — $(date)"
echo "  contigs   : ${ASM_DIR}/"
echo "  scaffolds : ${SCAFFOLD_DIR}/${SAMPLE}.hap{1,2}_scaffolds_final.fa"
echo "  QUAST     : ${QC_DIR}/quast/"
echo "  Juicebox  : ${SCAFFOLD_DIR}/${SAMPLE}.hap{1,2}_JBAT.{hic,assembly}"
echo ""
echo "  Next, once every sample has finished:"
echo "    bash   03_ptarmigan_reference.sh"
echo "    sbatch 04_chromosome_homology.sh"
echo "    N=\$(grep -v '^#' assembly_manifest.tsv | tail -n +2 | grep -c .)"
echo "    sbatch --array=0-\$((N-1))%6 05_ragtag_liftoff_array.sh"
echo "============================================================"
