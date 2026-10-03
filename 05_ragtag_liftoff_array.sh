#!/bin/bash
# =============================================================================
# SLURM ARRAY JOB: RE-SCAFFOLD AND RE-ANNOTATE AGAINST ROCK PTARMIGAN
# Step 05 — requires 03_ptarmigan_reference.sh, and requires
# 02_assembly_scaffold_array.sh to have produced the yahs scaffolds.
#
# This is the reference-guided half of the pipeline. hifiasm and yahs are NOT
# re-run: the Hi-C-scaffolded assemblies step 02 produced are reference-free,
# and are the correct input to a reference-guided ordering. Step 02's tree is
# read, never written to, so this step can be re-run against a different
# reference without rebuilding anything upstream of it.
#
# INPUT  (read-only, from step 02)
#   ${PROJECT_DIR}/yahs/<SAMPLE>.<HAP>_scaffolds_final.fa
#
# OUTPUT (a parallel tree, suffixed _ptarmigan)
#   ragtag_ptarmigan/<SAMPLE>.<HAP>/ragtag.scaffold.{fasta,agp}
#   ragtag_ptarmigan/<SAMPLE>.<HAP>/<SAMPLE>.<HAP>.rename_map.tsv
#   final_ptarmigan/<SPECIES>_<SAMPLE>_<HAP>.pseudo_chr.fasta
#   final_ptarmigan/<SPECIES>_<SAMPLE>_<HAP>.unplaced_short.fasta
#   final_ptarmigan/<SPECIES>_<SAMPLE>_<HAP>.liftoff.gff3
#   qc_ptarmigan/<SAMPLE>.<HAP>.hic2final.cram        (if DO_HIC=true)
#   qc_ptarmigan/<SAMPLE>.<HAP>.pretext               (if DO_HIC=true)
#   qc_ptarmigan/quast/<SAMPLE>.<HAP>.post_ragtag/
#   qc_ptarmigan/depth_check/<SAMPLE>.<HAP>.coverage.tsv   (if DO_DEPTH=true)
#
# WHY THE HI-C IS RE-ALIGNED (DO_HIC)
#   Hi-C reads are aligned against the FINAL, reference-ordered assembly here,
#   not reused from step 02. Step 02's CRAM is aligned to the pre-RagTag yahs
#   scaffolds; the pseudo-chromosome assembly is the same sequence in a
#   different order under different names, so those coordinates are meaningless
#   against it and reusing that CRAM would silently produce nonsense
#   join-support scores. Step 11 therefore needs its own alignment, which is
#   what DO_HIC produces. It is also the expensive part of this job: with
#   DO_HIC=false the job is a few hours, with it on, most of two days. Turn it
#   off if you only want the karyotype and BUSCO results; turn it on before
#   running 11_join_support.sh.
#
# USAGE
#   N=$(grep -v '^#' assembly_manifest.tsv | tail -n +2 | grep -c .)
#   sbatch --array=0-$((N-1))%6 05_ragtag_liftoff_array.sh
#
#   Rerun one failed task, e.g. index 7:
#   sbatch --array=7 05_ragtag_liftoff_array.sh
#
# Every stage is resume-safe and timestamp-checked (`-nt`) against its input,
# so a resubmitted task skips completed work and a rebuilt input regenerates
# everything downstream of it rather than keeping stale output.
# =============================================================================
#SBATCH --job-name=grouse_ragtag_ptarmigan
#SBATCH --output=logs/%x_%A_%a.out
#SBATCH --error=logs/%x_%A_%a.err
#SBATCH -A dewoody
# 3 days covers both haplotypes with DO_HIC=true AND DO_DEPTH=true: bwa mem over
# the Hi-C library dominates, and the HiFi realignment for the depth check adds a
# second full-library pass. With both off this finishes in well under a day —
# trim the walltime then, it will schedule sooner.
#SBATCH -t 3-00:00:00
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=24
# 96G: Liftoff is the memory peak (it holds the feature database and minimap2
# index together), not bwa or RagTag. Haplotypes run sequentially so the peak is
# per-haplotype. Check `sacct -o MaxRSS` after the first few and trim.
#SBATCH --mem=96G
#SBATCH -p cpu
#SBATCH --mail-type=BEGIN,END,FAIL
#SBATCH --mail-user=blackan@purdue.edu

# =============================================================================
# ENVIRONMENT SETUP
# =============================================================================
set -euo pipefail

# Defensive: --export=ALL (the sbatch default) carries the submitting shell's
# Lmod state into the job, and liftoff refuses to load alongside anaconda.
module unload anaconda 2>/dev/null || true

ml biocontainers
ml quast
ml bwa
ml samtools/1.22.1

# unset LD_PRELOAD before anything else runs: RCAC's XALT usage-tracking library
# is injected this way and fails on some nodes (GLIBC_2.33/2.34 mismatch), which
# can kill subshells under `set -e`. It's accounting only — safe to drop.
unset LD_PRELOAD

# Unlike step 02 these can be loaded up front: this job never touches anaconda
# (PretextMap and the dotplot interpreter are invoked by absolute path from the
# conda envs step 02 already built), so the Lmod conflict that forced step 02 to
# defer these does not arise here.
ml ragtag
ml liftoff
ml minimap2

# =============================================================================
# USER-DEFINED VARIABLES
# =============================================================================
MANIFEST="${SLURM_SUBMIT_DIR}/assembly_manifest.tsv"

PROJECT_DIR="${CLUSTER_SCRATCH}/GROUSE/grouse_asm"
REF_DIR="${PROJECT_DIR}/ref"

# ---- Read-only inputs from the existing tree --------------------------------
SCAFFOLD_DIR="${PROJECT_DIR}/yahs"          # yahs output — the real input here
FILT_DIR="${PROJECT_DIR}/hifi_filtered"     # adapter-filtered HiFi from step 02
CONDA_ENVS_DIR="${PROJECT_DIR}/conda_envs"  # envs built by step 02

# ---- Parallel output tree ---------------------------------------------------
# Changing REF_TAG alone retargets every output path, which is what makes a
# third reference (turkey, say) a one-line change rather than a new script.
REF_TAG="ptarmigan"
RAGTAG_DIR="${PROJECT_DIR}/ragtag_${REF_TAG}"
LIFTOFF_DIR="${PROJECT_DIR}/liftoff_${REF_TAG}"
FINAL_DIR="${PROJECT_DIR}/final_${REF_TAG}"
QC_DIR="${PROJECT_DIR}/qc_${REF_TAG}"

# ---- Reference --------------------------------------------------------------
# Resolved and recorded by step 03. Read rather than hardcoded so the two
# scripts cannot disagree about which assembly was downloaded.
ASM_DIR_RECORD="${REF_DIR}/ptarmigan_asm_dir.txt"

# ---- What chr_N is going to MEAN --------------------------------------------
# native   (DEFAULT) use the ptarmigan assembly report's own names. For bLagMut1
#            those are SUPER_N stripped to N. Be clear about what that is: a
#            LENGTH RANKING, not homology, so chr_6 here is NOT chicken chr6.
#            These are the names the delivered assemblies carry. They name a
#            chromosome rather than asserting a comparison, which is what you
#            want in a deliverable; the correspondence to the chicken karyotype
#            is recorded separately by step 04 and is the thing to cite when
#            relating these assemblies to the galliform literature.
# homology : use the chicken-homology map from 04_chromosome_homology.sh, so
#            chr_N names the same chromosome as chicken chr_N. A ptarmigan
#            chromosome spanning two chicken chromosomes is then named for both
#            (chr_6_8), one spanning part of one is lettered (chr_2a), and one
#            with no confident match becomes chr_u<N>. Useful while comparing
#            references; poor as an identifier, because every name encodes a
#            claim about a third genome. REQUIRES step 04 to have run.
#            06_rename_chromosomes.sh converts a tree built this way to native
#            names after the fact, which is how the delivered assemblies were
#            produced.
#
# Overridable at submit time: --export=ALL,CHR_NAMING=homology
CHR_NAMING="${CHR_NAMING:-native}"

HOMOLOGY_MAP="${PROJECT_DIR}/qc/chromosome_homology/ptarmigan_vs_chicken.rename_map.tsv"

# ---- Stage toggles ----------------------------------------------------------
# Each can also be set at submit time without editing this file, e.g.
#   sbatch --export=ALL,DO_HIC=false --array=0-22%6 05_ragtag_liftoff_array.sh
# which is the convenient way to do a quick RagTag-only pass first and add the
# expensive Hi-C stage later.
DO_LIFTOFF="${DO_LIFTOFF:-true}"
DO_QUAST="${DO_QUAST:-true}"

# Hi-C re-alignment + Pretext contact map. Required by step 11; see the header.
DO_HIC="${DO_HIC:-true}"

# Orientation dotplot vs the ptarmigan reference (independent minimap2 pass,
# i.e. not RagTag's own alignment). Cheap and it is the fastest visual check
# that the new ordering is sane.
DO_DOTPLOT="${DO_DOTPLOT:-true}"

# tidk telomere profiling is OFF by default here, unlike step 02. It measures
# telomeric repeat density along each chromosome, which depends on sequence
# content rather than on the reference used to order it — so re-running it on
# the same scaffolds in a different order reproduces the step 02 result almost
# exactly. Set true if you want the figure regenerated against the new names.
DO_TIDK="${DO_TIDK:-false}"

# HiFi depth check: realign the adapter-filtered HiFi reads to the FINAL
# assembly and compare chr_Z / chr_W depth against the autosome mean. This is
# the read-level corroboration of the sex calls that 10_sex_chromosome_check.sh
# makes from assembly content and BUSCO placement: a female haplotype carrying
# a real W shows chr_W near 1x of the per-haplotype autosome mean, while a
# haplotype that simply failed to assemble its Z shows chr_Z far below it. Two
# independent lines of evidence for the same claim is the point. Costs a
# minimap2 pass over the full HiFi set per haplotype — hours, not minutes — so
# turn it off for a quick structural-only pass.
DO_DEPTH="${DO_DEPTH:-true}"

HIC_MIN_MAPQ=20

# chr_* sequences are
# always kept regardless of length — galliform microchromosomes are legitimately
# small — and short unplaced scaffolds are set aside, not discarded.
MIN_UNPLACED_SCAFFOLD_LEN=50000

# Absolute-path tools from step 02's conda envs.
PRETEXT_ENV_DIR="${CONDA_ENVS_DIR}/pretext"
PRETEXTMAP_BIN="${PRETEXT_ENV_DIR}/bin/PretextMap"
PRETEXTSNAPSHOT_BIN="${PRETEXT_ENV_DIR}/bin/PretextSnapshot"
DOTPLOT_PYTHON_BIN="${CONDA_ENVS_DIR}/dotplot-python/bin/python3"
TIDK_BIN="${CONDA_ENVS_DIR}/tidk/bin/tidk"

THREADS=$SLURM_CPUS_PER_TASK

mkdir -p logs "$RAGTAG_DIR" "$LIFTOFF_DIR" "$FINAL_DIR" "$QC_DIR" \
    "${QC_DIR}/quast" "${QC_DIR}/tidk"

# =============================================================================
# PRE-FLIGHT
# =============================================================================
for BIN in ragtag.py liftoff minimap2 bwa samtools quast.py; do
    if ! command -v "$BIN" >/dev/null 2>&1; then
        echo "ERROR: '${BIN}' not on PATH after module load."
        echo "       Check: ml spider ${BIN%%.*}"
        exit 1
    fi
done

if [[ ! -s "$ASM_DIR_RECORD" ]]; then
    echo "ERROR: ${ASM_DIR_RECORD} not found."
    echo "       Run 03_ptarmigan_reference.sh first."
    exit 1
fi
ASM_DIR="$(head -n1 "$ASM_DIR_RECORD")"
if [[ -z "$ASM_DIR" ]]; then
    echo "ERROR: ${ASM_DIR_RECORD} is empty."
    exit 1
fi

REF_FASTA="${REF_DIR}/${ASM_DIR}_genomic.fna"
REF_GFF="${REF_DIR}/${ASM_DIR}_genomic.gff"
REF_CHR_MAP="${REF_DIR}/${ASM_DIR}.chr_map.tsv"

echo ">>> Reference : ${ASM_DIR}"
echo ">>> FASTA     : ${REF_FASTA}"

if [[ ! -s "$REF_FASTA" ]]; then
    echo "ERROR: reference FASTA not found: ${REF_FASTA}"
    echo "       Run 03_ptarmigan_reference.sh first."
    exit 1
fi
if [[ ! -s "$REF_CHR_MAP" ]]; then
    echo "ERROR: chromosome map not found: ${REF_CHR_MAP}"
    echo "       Run 03_ptarmigan_reference.sh first."
    exit 1
fi

# Choose the map that decides what chr_N means. This is the single most
# consequential setting in the script, so it is resolved once, loudly, up front
# rather than inside the per-haplotype loop.
case "$CHR_NAMING" in
    homology)
        if [[ ! -s "$HOMOLOGY_MAP" ]]; then
            echo "ERROR: CHR_NAMING=homology but the homology map is missing:"
            echo "         ${HOMOLOGY_MAP}"
            echo ""
            echo "       Run 04_chromosome_homology.sh first. It aligns the"
            echo "       ptarmigan reference to chicken and works out which"
            echo "       chicken chromosome each ptarmigan chromosome corresponds"
            echo "       to."
            echo ""
            echo "       To use the reference's own numbering instead (the"
            echo "       default for this pipeline):"
            echo "         sbatch --export=ALL,CHR_NAMING=native ..."
            exit 1
        fi
        ACTIVE_CHR_MAP="$HOMOLOGY_MAP"
        echo ">>> chr naming: HOMOLOGY — chr_N means chicken chr_N"
        echo ">>> chr map   : ${HOMOLOGY_MAP}"
        echo "                $(grep -vc '^#' "$HOMOLOGY_MAP") chromosomes"
        if grep -qv '^#' <<< "$(grep '_' "$HOMOLOGY_MAP" | grep -v '^#' || true)"; then
            echo ">>> compound names present (a ptarmigan chromosome spanning"
            echo "    several chicken chromosomes), which is the expected shape"
            echo "    if the karyotypes genuinely differ:"
            awk -F'\t' '!/^#/ && $2 ~ /_/ {printf "      chr_%s\n", $2}' "$HOMOLOGY_MAP" | head -6
        fi
        ;;
    native)
        ACTIVE_CHR_MAP="$REF_CHR_MAP"
        echo ">>> chr naming: NATIVE — the ptarmigan assembly report's own names"
        echo ">>> chr map   : $(wc -l < "$REF_CHR_MAP") assembled molecules"
        echo "    !! For bLagMut1 these are a LENGTH RANKING, not homology."
        echo "       chr_N here is NOT chicken chr_N. To relate these assemblies"
        echo "       to the galliform literature use the correspondence table"
        echo "       from 04_chromosome_homology.sh, never the number alone."
        ;;
    *)
        echo "ERROR: CHR_NAMING must be 'homology' or 'native', got '${CHR_NAMING}'"
        exit 1
        ;;
esac

if [[ "$DO_LIFTOFF" == true && ! -s "$REF_GFF" ]]; then
    echo ""
    echo "  !! DO_LIFTOFF=true but no annotation GFF at:"
    echo "     ${REF_GFF}"
    echo "     Step 03 reported whether one was available. Continuing with"
    echo "     Liftoff disabled — RagTag, BUSCO, karyotype and join support are"
    echo "     all unaffected."
    DO_LIFTOFF=false
fi

if [[ "$DO_HIC" == true ]]; then
    for BIN in "$PRETEXTMAP_BIN" "$PRETEXTSNAPSHOT_BIN"; do
        if [[ ! -x "$BIN" ]]; then
            echo "ERROR: DO_HIC=true but ${BIN} is missing."
            echo "       It is created by 02_assembly_scaffold_array.sh's one-time"
            echo "       setup. Either run that once, or set DO_HIC=false."
            exit 1
        fi
    done
fi
if [[ "$DO_DOTPLOT" == true && ! -x "$DOTPLOT_PYTHON_BIN" ]]; then
    echo "  !! DOTPLOT_PYTHON_BIN missing (${DOTPLOT_PYTHON_BIN}) — skipping dotplots"
    DO_DOTPLOT=false
fi
if [[ "$DO_TIDK" == true && ! -x "$TIDK_BIN" ]]; then
    echo "  !! TIDK_BIN missing (${TIDK_BIN}) — skipping tidk"
    DO_TIDK=false
fi

# =============================================================================
# RESOLVE SAMPLE FOR THIS ARRAY TASK
# =============================================================================
if [[ -z "${SLURM_ARRAY_TASK_ID:-}" ]]; then
    echo "ERROR: SLURM_ARRAY_TASK_ID is not set."
    echo "Submit with: sbatch --array=0-N 05_ragtag_liftoff_array.sh"
    exit 1
fi
if [[ ! -f "$MANIFEST" ]]; then
    echo "ERROR: Manifest not found: ${MANIFEST}"
    exit 1
fi

mapfile -t ROWS < <(grep -v '^#' "$MANIFEST" | tail -n +2 | grep -v '^[[:space:]]*$')
LINE="${ROWS[$SLURM_ARRAY_TASK_ID]:-}"
if [[ -z "$LINE" ]]; then
    echo "ERROR: No manifest row at index ${SLURM_ARRAY_TASK_ID} (${#ROWS[@]} samples)"
    exit 1
fi

# Same column layout as step 02: sample, species, hifi bams, hic R1, hic R2, ONT.
IFS=$'\t' read -r SAMPLE SPECIES _HIFI_BAMS HIC_R1 HIC_R2 _ONT_UL <<< "$LINE"
if [[ -z "${SAMPLE:-}" || -z "${SPECIES:-}" ]]; then
    echo "ERROR: Malformed manifest row at index ${SLURM_ARRAY_TASK_ID}: ${LINE}"
    exit 1
fi

if [[ "$DO_HIC" == true ]]; then
    if [[ -z "${HIC_R1:-}" || -z "${HIC_R2:-}" ]]; then
        echo "ERROR: DO_HIC=true but the manifest row has no Hi-C reads:"
        echo "       ${LINE}"
        exit 1
    fi
    for F in "$HIC_R1" "$HIC_R2"; do
        if [[ ! -s "$F" ]]; then
            echo "ERROR: Hi-C read file not found: ${F}"
            exit 1
        fi
    done
fi

echo ""
echo ">>> Array task ${SLURM_ARRAY_TASK_ID} -> ${SAMPLE} (${SPECIES})"
echo ">>> Started: $(date)"
echo ">>> Threads: ${THREADS}"
echo ">>> Stages : ragtag=yes liftoff=${DO_LIFTOFF} quast=${DO_QUAST} hic=${DO_HIC} dotplot=${DO_DOTPLOT} tidk=${DO_TIDK}"

# =============================================================================
# PER-HAPLOTYPE
# =============================================================================
for HAP in hap1 hap2; do

    PREFIX="${SPECIES}_${SAMPLE}_${HAP}"
    SCAFFOLDS="${SCAFFOLD_DIR}/${SAMPLE}.${HAP}_scaffolds_final.fa"

    echo ""
    echo "============================================================"
    echo ">>> ${PREFIX}"
    echo "============================================================"

    if [[ ! -s "$SCAFFOLDS" ]]; then
        echo "ERROR: yahs scaffolds not found: ${SCAFFOLDS}"
        echo "       02_assembly_scaffold_array.sh must have completed for this"
        echo "       sample/haplotype. This script re-scaffolds its output; it"
        echo "       does not re-run hifiasm or yahs."
        exit 1
    fi

    # ---------------------------------------------------------------- Step 1
    echo ""
    echo ">>> [1] RagTag vs ${ASM_DIR}"
    RAGTAG_OUT="${RAGTAG_DIR}/${SAMPLE}.${HAP}"
    PSEUDO_CHR="${RAGTAG_OUT}/ragtag.scaffold.fasta"

    if [[ -s "$PSEUDO_CHR" && "$PSEUDO_CHR" -nt "$SCAFFOLDS" ]]; then
        echo "  RagTag output already exists and is current — skipping"
    else
        # -u writes the unplaced sequences into the output as well, matching
        # step 02 exactly: nothing yahs produced is dropped at this stage.
        ragtag.py scaffold -o "$RAGTAG_OUT" -t "$THREADS" -u "$REF_FASTA" "$SCAFFOLDS"
    fi
    if [[ ! -s "$PSEUDO_CHR" ]]; then
        echo "ERROR: RagTag did not produce expected output: ${PSEUDO_CHR}"
        exit 1
    fi
    if [[ ! -s "${RAGTAG_OUT}/ragtag.scaffold.agp" ]]; then
        echo "ERROR: RagTag produced no AGP: ${RAGTAG_OUT}/ragtag.scaffold.agp"
        echo "       Step 11 needs it to score the joins."
        exit 1
    fi

    # ---------------------------------------------------------------- Step 2
    echo ""
    echo ">>> [2] Rename to ptarmigan chromosome names, strip _RagTag"
    FINAL_FASTA="${FINAL_DIR}/${PREFIX}.pseudo_chr.fasta"
    SHORT_FASTA="${FINAL_DIR}/${PREFIX}.unplaced_short.fasta"
    RENAME_MAP="${RAGTAG_OUT}/${SAMPLE}.${HAP}.rename_map.tsv"

    # Guard requires BOTH that renaming happened AND that the length filter ran
    # (the .unplaced_short.fasta companion exists), and that the result is newer
    # than the RagTag output it came from.
    ALREADY_RENAMED=false
    if [[ -s "$FINAL_FASTA" && -f "$SHORT_FASTA" && "$FINAL_FASTA" -nt "$PSEUDO_CHR" ]] \
        && grep -qm1 '^>chr_\|^>scaffold_' "$FINAL_FASTA" \
        && ! grep -qm1 '_RagTag' "$FINAL_FASTA"; then
        ALREADY_RENAMED=true
    fi

    if [[ "$ALREADY_RENAMED" == true ]]; then
        echo "  already renamed — skipping"
    else
        samtools faidx "$PSEUDO_CHR"

        # Names that pass through the rewrite unchanged: RagTag leaves UNPLACED
        # input sequences under their original yahs names (scaffold_1,
        # scaffold_2, ...), and only sequences it placed get a _RagTag suffix.
        # The generated scaffold_N numbering below must avoid these, or two
        # different sequences end up with the same name — see the comment there.
        TAKEN_NAMES="${RAGTAG_OUT}/${SAMPLE}.${HAP}.passthrough_names.txt"
        awk -F'\t' '$1 !~ /_RagTag$/ {print $1}' "${PSEUDO_CHR}.fai" > "$TAKEN_NAMES"

        # Placed sequences carry a _RagTag suffix and are named for the
        # reference accession they were assigned to; map those to chr_<name>
        # via the assembly report. Everything else becomes scaffold_N numbered
        # by descending length, so the numbering is stable across reruns.
        #
        # The `do ... while` skips any scaffold_N that a pass-through sequence
        # already owns. Without it, a grouse scaffold that RagTag placed onto one
        # of the REFERENCE's unplaced scaffolds (an NW_ accession, absent from
        # the chromosome map) is named scaffold_1, while a yahs scaffold that
        # RagTag could not place at all is also called scaffold_1 — two distinct
        # sequences, one name. samtools faidx rejects the duplicate, so it
        # surfaces as an opaque indexing failure rather than as the naming
        # conflict it is. Whether it bites depends entirely on whether RagTag
        # placed anything onto the reference's unplaced scaffolds, so a reference
        # with a larger unplaced fraction makes it likely rather than impossible.
        # When there is no clash this produces exactly the same names as before.
        awk -F'\t' '$1 ~ /_RagTag$/ {print $1"\t"$2}' "${PSEUDO_CHR}.fai" \
            | sort -k2,2 -nr \
            | awk -F'\t' -v chrmap="$ACTIVE_CHR_MAP" -v takenfile="$TAKEN_NAMES" '
                BEGIN {
                    # Skip comments and blanks: the homology map from step 04
                    # carries a "# accession  chromosome_name" header, and the
                    # plain assembly-report map does not. Both are read here.
                    while ((getline line < chrmap) > 0) {
                        if (line ~ /^#/ || line ~ /^[[:space:]]*$/) continue
                        split(line, a, "\t")
                        if (a[1] != "" && a[2] != "") chrname[a[1]] = a[2]
                    }
                    while ((getline line < takenfile) > 0) {
                        taken[line] = 1
                    }
                }
                {
                    acc = $1
                    sub(/_RagTag$/, "", acc)
                    if (acc in chrname) {
                        print $1"\tchr_"chrname[acc]
                    } else {
                        do { scafn++ } while (("scaffold_" scafn) in taken)
                        print $1"\tscaffold_"scafn
                    }
                }
            ' > "$RENAME_MAP"

        if [[ ! -s "$RENAME_MAP" ]]; then
            echo "ERROR: rename map came out empty for ${PREFIX}"
            echo "       No _RagTag sequences in ${PSEUDO_CHR}.fai?"
            exit 1
        fi

        # If nothing mapped to a real chromosome, the chr map and the RagTag
        # output are keyed on different accession styles and every sequence has
        # silently become scaffold_N. That is a wrong-reference bug, not a
        # result, so stop here rather than producing 46 unusable assemblies.
        N_CHR_ASSIGNED=$(awk -F'\t' '$2 ~ /^chr_/' "$RENAME_MAP" | wc -l)
        if (( N_CHR_ASSIGNED == 0 )); then
            echo "ERROR: no sequence was assigned a chromosome name for ${PREFIX}."
            echo "       The RagTag sequence names and ${ACTIVE_CHR_MAP} do not"
            echo "       share an accession style. First few of each:"
            echo "       --- ragtag.scaffold.fasta.fai ---"
            head -3 "${PSEUDO_CHR}.fai" | cut -f1 | sed 's/^/         /'
            echo "       --- chr_map.tsv ---"
            head -3 "$ACTIVE_CHR_MAP" | sed 's/^/         /'
            exit 1
        fi
        echo "  rename map: ${RENAME_MAP} (${N_CHR_ASSIGNED} chromosomes)"

        awk -v map="$RENAME_MAP" '
            BEGIN {
                while ((getline line < map) > 0) {
                    split(line, a, "\t")
                    newname[a[1]] = a[2]
                }
            }
            /^>/ {
                split(substr($0, 2), parts, " ")
                name = parts[1]
                if (name in newname) {
                    print ">" newname[name]
                } else {
                    sub(/_RagTag$/, "", name)
                    print ">" name
                }
                next
            }
            { print }
        ' "$PSEUDO_CHR" > "${FINAL_FASTA}.prefilter"

        # Independent check that no two sequences ended up sharing a name. The
        # skip-taken-numbers logic above should make this impossible; this
        # verifies it rather than trusting it, because the failure mode is a
        # silently corrupt assembly and the cost of checking is one sort.
        DUP_NAMES=$(grep '^>' "${FINAL_FASTA}.prefilter" | sed 's/^>//' \
                    | awk '{print $1}' | sort | uniq -d || true)
        if [[ -n "$DUP_NAMES" ]]; then
            echo "ERROR: duplicate sequence names after renaming ${PREFIX}:"
            printf '%s\n' "$DUP_NAMES" | sed 's/^/         /'
            echo "       This should not be reachable. Inspect ${RENAME_MAP}"
            echo "       and ${TAKEN_NAMES}."
            exit 1
        fi

        # --- Minimum-length filter on UNPLACED scaffolds only ---
        samtools faidx "${FINAL_FASTA}.prefilter"
        awk -F'\t' -v min="$MIN_UNPLACED_SCAFFOLD_LEN" \
            '$1 ~ /^chr_/ || $2 >= min {print $1}' "${FINAL_FASTA}.prefilter.fai" > "${FINAL_FASTA}.keep.txt"
        awk -F'\t' -v min="$MIN_UNPLACED_SCAFFOLD_LEN" \
            '$1 !~ /^chr_/ && $2 < min {print $1}' "${FINAL_FASTA}.prefilter.fai" > "${FINAL_FASTA}.short.txt"

        if [[ ! -s "${FINAL_FASTA}.keep.txt" ]]; then
            echo "ERROR: length filter would keep nothing — check ${FINAL_FASTA}.prefilter"
            exit 1
        fi
        samtools faidx "${FINAL_FASTA}.prefilter" -r "${FINAL_FASTA}.keep.txt" > "$FINAL_FASTA"
        if [[ -s "${FINAL_FASTA}.short.txt" ]]; then
            samtools faidx "${FINAL_FASTA}.prefilter" -r "${FINAL_FASTA}.short.txt" > "$SHORT_FASTA"
        else
            : > "$SHORT_FASTA"
        fi

        N_IN=$(wc -l < "${FINAL_FASTA}.prefilter.fai")
        N_KEEP=$(wc -l < "${FINAL_FASTA}.keep.txt")
        N_SHORT=$(wc -l < "${FINAL_FASTA}.short.txt")
        if (( N_KEEP + N_SHORT != N_IN )); then
            echo "ERROR: sequence counts do not reconcile for ${PREFIX}:"
            echo "       in ${N_IN}, kept ${N_KEEP}, short ${N_SHORT}"
            exit 1
        fi
        echo "  kept     : ${N_KEEP} sequences (all chr_* plus unplaced >= ${MIN_UNPLACED_SCAFFOLD_LEN} bp)"
        echo "  set aside: ${N_SHORT} short unplaced scaffolds -> $(basename "$SHORT_FASTA")"

        rm -f "${FINAL_FASTA}.prefilter" "${FINAL_FASTA}.prefilter.fai" \
              "${FINAL_FASTA}.keep.txt" "${FINAL_FASTA}.short.txt"

        # Liftoff caches a minimap2 index (<target>.mmi) by PATH, not content,
        # and bwa's index is equally stale once FINAL_FASTA is rewritten. Clear
        # both so nothing downstream reads an index for a different sequence.
        rm -f "${FINAL_FASTA}.mmi" "${FINAL_FASTA}".{amb,ann,bwt,pac,sa}
        rm -rf "${LIFTOFF_DIR}/${SAMPLE}.${HAP}_intermediate"
    fi

    samtools faidx "$FINAL_FASTA"
    echo "  Final assembly: ${FINAL_FASTA}"
    echo "  Chromosomes:"
    awk -F'\t' '$1 ~ /^chr_/ {printf "    %-10s %12.2f Mb\n", $1, $2/1e6}' "${FINAL_FASTA}.fai" | head -50

    # ---------------------------------------------------------------- Step 3
    if [[ "$DO_QUAST" == true ]]; then
        echo ""
        echo ">>> [3] QUAST — final assembly"
        QUAST_OUT="${QC_DIR}/quast/${SAMPLE}.${HAP}.post_ragtag"
        if [[ -f "${QUAST_OUT}/report.txt" && "${QUAST_OUT}/report.txt" -nt "$FINAL_FASTA" ]]; then
            echo "  already exists and is current — skipping"
        else
            quast.py -o "$QUAST_OUT" -t "$THREADS" --large "$FINAL_FASTA"
        fi
    fi

    # ---------------------------------------------------------------- Step 4
    if [[ "$DO_LIFTOFF" == true ]]; then
        echo ""
        echo ">>> [4] Liftoff ptarmigan annotation"
        LIFTOFF_GFF="${FINAL_DIR}/${PREFIX}.liftoff.gff3"
        LIFTOFF_UNMAPPED="${LIFTOFF_DIR}/${SAMPLE}.${HAP}.unmapped_features.txt"
        LIFTOFF_INTERMEDIATE="${LIFTOFF_DIR}/${SAMPLE}.${HAP}_intermediate"

        if [[ -s "$LIFTOFF_GFF" && "$LIFTOFF_GFF" -nt "$FINAL_FASTA" ]]; then
            echo "  already exists and is current — skipping"
        else
            mkdir -p "$LIFTOFF_INTERMEDIATE"
            liftoff \
                -g "$REF_GFF" \
                -o "$LIFTOFF_GFF" \
                -u "$LIFTOFF_UNMAPPED" \
                -dir "$LIFTOFF_INTERMEDIATE" \
                -p "$THREADS" \
                "$FINAL_FASTA" \
                "$REF_FASTA"
        fi
        echo "  Liftoff annotation: ${LIFTOFF_GFF}"
        if [[ -s "$LIFTOFF_UNMAPPED" ]]; then
            echo "  Unmapped features : $(wc -l < "$LIFTOFF_UNMAPPED") -> ${LIFTOFF_UNMAPPED}"
            echo "  (a closer reference transfers more features; this count is"
            echo "   worth reporting, and worth comparing against any other"
            echo "   reference you run)"
        fi
    fi

    # ---------------------------------------------------------------- Step 5
    if [[ "$DO_DOTPLOT" == true ]]; then
        echo ""
        echo ">>> [5] Orientation dotplot vs ptarmigan (independent minimap2)"
        PAF="${QC_DIR}/${SAMPLE}.${HAP}.vs_ptarmigan.paf"
        DOTPLOT="${QC_DIR}/${SAMPLE}.${HAP}.orientation_dotplot.png"

        if [[ -s "$DOTPLOT" && "$DOTPLOT" -nt "$FINAL_FASTA" ]]; then
            echo "  already exists and is current — skipping"
        else
            # asm20 as in step 02: these are between-genus comparisons, so the
            # divergence preset has to be permissive. Lagopus-Tympanuchus is
            # closer than Gallus-Tympanuchus, so asm20 is, if anything,
            # conservative here — keeping it identical is what makes the two
            # dotplots comparable.
            minimap2 -x asm20 -t "$THREADS" "$REF_FASTA" "$FINAL_FASTA" > "$PAF"
            "$DOTPLOT_PYTHON_BIN" - "$PAF" "$DOTPLOT" \
                "${SPECIES} ${SAMPLE} ${HAP} vs ptarmigan (${ASM_DIR})" << 'PYEOF'
import math, sys
from collections import defaultdict
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

paf_path, out_path, title = sys.argv[1:4]

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
        if match_len < 1000:
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

fig.suptitle(title, fontsize=10)
fig.tight_layout(rect=(0, 0, 1, 0.97))
fig.savefig(out_path, dpi=150)
print(f"Wrote {out_path}")
PYEOF
        fi
        echo "  Orientation dotplot: ${DOTPLOT}"
    fi

    # ---------------------------------------------------------------- Step 6
    if [[ "$DO_HIC" == true ]]; then
        echo ""
        echo ">>> [6] Hi-C re-alignment + Pretext contact map"
        FINAL_HIC_CRAM="${QC_DIR}/${SAMPLE}.${HAP}.hic2final.cram"
        PRETEXT_MAP="${QC_DIR}/${SAMPLE}.${HAP}.pretext"

        if [[ -s "$PRETEXT_MAP" && -s "$FINAL_HIC_CRAM" && "$PRETEXT_MAP" -nt "$FINAL_FASTA" ]]; then
            echo "  already exists and is current — skipping"
        else
            if [[ ! -f "${FINAL_FASTA}.bwt" || "$FINAL_FASTA" -nt "${FINAL_FASTA}.bwt" ]]; then
                echo "  bwa index"
                bwa index "$FINAL_FASTA"
            fi
            # Identical flags to step 02 so the two CRAMs — and therefore the
            # two sets of join-support scores — are produced the same way:
            # -5SP for Hi-C, MAPQ filter at HIC_MIN_MAPQ. CRAM is fine, it is
            # only ever read back through samtools view.
            echo "  bwa mem -5SP (this is the long step)"
            bwa mem -5SP -t "$THREADS" "$FINAL_FASTA" "$HIC_R1" "$HIC_R2" \
                | samtools view -@ "$THREADS" -buS -q "$HIC_MIN_MAPQ" - \
                | samtools sort -@ "$THREADS" -m 1G --output-fmt cram \
                    --reference "$FINAL_FASTA" -o "${FINAL_HIC_CRAM}.part" -
            mv "${FINAL_HIC_CRAM}.part" "$FINAL_HIC_CRAM"
            samtools index "$FINAL_HIC_CRAM"

            samtools view -h -T "$FINAL_FASTA" "$FINAL_HIC_CRAM" \
                | "$PRETEXTMAP_BIN" -o "${PRETEXT_MAP}.part" \
                    --sortby length --sortorder descend --mapq "$HIC_MIN_MAPQ"
            mv "${PRETEXT_MAP}.part" "$PRETEXT_MAP"

            # --printSequenceNames labels each sequence on the map, which is
            # what makes a 40-chromosome contact map readable at all. Older
            # PretextSnapshot builds lack the flag, hence the fallback below.
            "$PRETEXTSNAPSHOT_BIN" --map "$PRETEXT_MAP" --sequences "=full" \
                --printSequenceNames --prefix "${SAMPLE}.${HAP}." --folder "$QC_DIR" \
                || "$PRETEXTSNAPSHOT_BIN" --map "$PRETEXT_MAP" --sequences "=full" \
                    --prefix "${SAMPLE}.${HAP}." --folder "$QC_DIR"
        fi
        echo "  Hi-C CRAM       : ${FINAL_HIC_CRAM}"
        echo "  Contact map     : ${QC_DIR}/${SAMPLE}.${HAP}.*.png"
    fi

    # ---------------------------------------------------------------- Step 7
    if [[ "$DO_TIDK" == true ]]; then
        echo ""
        echo ">>> [7] tidk telomere profiling (named chromosomes only)"
        TIDK_PREFIX="${SAMPLE}.${HAP}"
        TIDK_CHR_FASTA="${QC_DIR}/tidk/${TIDK_PREFIX}.chr_only.fasta"
        TIDK_TSV="${QC_DIR}/tidk/${TIDK_PREFIX}_telomeric_repeat_windows.tsv"
        TIDK_PLOT="${QC_DIR}/tidk/${TIDK_PREFIX}.svg"

        if [[ -s "$TIDK_PLOT" && "$TIDK_PLOT" -nt "$FINAL_FASTA" ]]; then
            echo "  already exists and is current — skipping"
        else
            awk -F'\t' '$1 ~ /^chr_/ {print $1}' "${FINAL_FASTA}.fai" > "${TIDK_CHR_FASTA}.names"
            if [[ ! -s "${TIDK_CHR_FASTA}.names" ]]; then
                echo "  WARNING: no chr_ sequences — skipping tidk"
                rm -f "${TIDK_CHR_FASTA}.names"
            else
                samtools faidx "$FINAL_FASTA" -r "${TIDK_CHR_FASTA}.names" > "$TIDK_CHR_FASTA"
                rm -f "${TIDK_CHR_FASTA}.names"
                "$TIDK_BIN" search --string TTAGGG --output "$TIDK_PREFIX" \
                    --dir "${QC_DIR}/tidk" "$TIDK_CHR_FASTA"
                "$TIDK_BIN" plot --tsv "$TIDK_TSV" --output "${QC_DIR}/tidk/${TIDK_PREFIX}"
                echo "  Telomere plot: ${TIDK_PLOT}"
            fi
        fi
    fi

    # ---------------------------------------------------------------- Step 8
    if [[ "$DO_DEPTH" == true ]]; then
        echo ""
        echo ">>> [8] HiFi depth check — chr_Z/chr_W vs autosomes"

        # Step 02 wrote these. Collect them without nullglob side effects and
        # without letting an empty glob reach minimap2 as a literal path.
        DEPTH_READS=()
        while IFS= read -r _r; do
            DEPTH_READS+=("$_r")
        done < <(find "${FILT_DIR}/${SAMPLE}" -maxdepth 1 -name '*.filt.fastq.gz' \
                      -type f 2>/dev/null | sort)

        if [[ ${#DEPTH_READS[@]} -eq 0 ]]; then
            echo "  WARNING: no adapter-filtered HiFi reads under ${FILT_DIR}/${SAMPLE}"
            echo "           Step 02 produces these. Skipping the depth check;"
            echo "           10_sex_chromosome_check.sh still runs without it."
        else
            mkdir -p "${QC_DIR}/depth_check"
            DEPTH_CRAM="${QC_DIR}/depth_check/${SAMPLE}.${HAP}.hifi2final.cram"
            DEPTH_TSV="${QC_DIR}/depth_check/${SAMPLE}.${HAP}.coverage.tsv"
            echo "  Reads: ${#DEPTH_READS[@]} file(s)"

            if [[ -s "$DEPTH_CRAM" && "$DEPTH_CRAM" -nt "$FINAL_FASTA" ]]; then
                echo "  HiFi-to-final alignment already current — skipping realignment"
            else
                minimap2 -ax map-hifi -t "$THREADS" "$FINAL_FASTA" "${DEPTH_READS[@]}" \
                    | samtools sort -@ "$THREADS" -m 1G \
                        --reference "$FINAL_FASTA" -O cram \
                        -o "${DEPTH_CRAM}.part" -
                mv "${DEPTH_CRAM}.part" "$DEPTH_CRAM"
                samtools index "$DEPTH_CRAM"
            fi

            samtools coverage --reference "$FINAL_FASTA" "$DEPTH_CRAM" > "${DEPTH_TSV}.part"
            mv "${DEPTH_TSV}.part" "$DEPTH_TSV"
            echo "  Per-sequence coverage: ${DEPTH_TSV}"

            # Baseline = mean depth over every named non-sex, non-organelle
            # chromosome. Matching on the chr_ prefix rather than chr_<digits>
            # keeps compound and lettered names (chr_6_8, chr_2a, chr_u30) in the
            # baseline when CHR_NAMING=homology.
            #
            # This is PER-HAPLOTYPE depth, so a genuine single-copy sequence —
            # a female's W in the haplotype that carries it — is expected near
            # 1.0x of the baseline, NOT 0.5x. The diagnostic signal is a value
            # near 0, meaning the sequence is barely present at all.
            awk -F'\t' '
                NR == 1 { next }
                $1 !~ /^chr_/ { next }
                $1 == "chr_Z" { z = $7; next }
                $1 == "chr_W" { w = $7; next }
                $1 == "chr_MT" { mt = $7; next }
                { auto_sum += $7; auto_n++ }
                END {
                    if (auto_n == 0) {
                        print "  no autosomes found to baseline against"
                        exit
                    }
                    m = auto_sum / auto_n
                    printf "  autosome mean depth (n=%d): %.2fx\n", auto_n, m
                    if (z != "") printf "  chr_Z : %7.2fx   ratio %.2f\n", z, z / m
                    if (w != "") printf "  chr_W : %7.2fx   ratio %.2f\n", w, w / m
                    if (mt != "") printf "  chr_MT: %7.2fx   ratio %.2f\n", mt, mt / m
                }
            ' "$DEPTH_TSV"
        fi
    fi

done

# =============================================================================
# DONE
# =============================================================================
echo ""
echo "============================================================"
echo ">>> ${SAMPLE} (${SPECIES}) complete: $(date)"
echo "  assemblies : ${FINAL_DIR}/${SPECIES}_${SAMPLE}_hap{1,2}.pseudo_chr.fasta"
if [[ "$DO_LIFTOFF" == true ]]; then
    echo "  annotation : ${FINAL_DIR}/${SPECIES}_${SAMPLE}_hap{1,2}.liftoff.gff3"
fi
echo "  AGP + map  : ${RAGTAG_DIR}/${SAMPLE}.hap{1,2}/"
if [[ "$DO_HIC" == true ]]; then
    echo "  Hi-C CRAM  : ${QC_DIR}/${SAMPLE}.hap{1,2}.hic2final.cram"
fi
if [[ "$DO_DEPTH" == true ]]; then
    echo "  Depth      : ${QC_DIR}/depth_check/${SAMPLE}.hap{1,2}.coverage.tsv"
fi
echo ""
echo "  Next, once every task has finished:"
echo "    N=\$(grep -v '^#' assembly_manifest.tsv | tail -n +2 | grep -c .)"
if [[ "$CHR_NAMING" == "homology" ]]; then
    echo "    sbatch 06_rename_chromosomes.sh   # homology names -> native names"
fi
echo "    sbatch --array=0-\$((N-1))%8 07_busco_array.sh"
echo "    sbatch 08_ptarmigan_reference_busco.sh"
echo "    sbatch 09_busco_plots.sh          # after 07 and 08"
echo "    sbatch 10_sex_chromosome_check.sh # after 07 and 08"
if [[ "$DO_HIC" == true ]]; then
    echo "    sbatch --array=0-\$((N-1))%8 11_join_support.sh"
fi
echo "============================================================"
