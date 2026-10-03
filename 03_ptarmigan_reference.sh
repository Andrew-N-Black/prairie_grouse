#!/bin/bash
# =============================================================================
# DOWNLOAD THE ROCK PTARMIGAN REFERENCE (GCF_023343835.1) + ANNOTATION
# Step 03 — the reference this pipeline actually assembles against.
#
# WHY A SECOND REFERENCE
#   Every assembly in this project was ordered against chicken (GRCg7b), a
#   Phasianinae genome roughly 37 My diverged from Tympanuchus. Two findings so
#   far are statements about chicken as much as about grouse:
#     - chicken chr6 + chr8 are one chromosome in all 46 haplotypes
#     - the RagTag join at 19.1-19.5 Mb on chr_4 (chicken GGA4p) is
#       Hi-C-unsupported in 28 assemblies across all three species
#   Both are consistent with the grouse retaining the ancestral galliform
#   karyotype while chicken carries derived fusions. The way to separate "our
#   assemblies are wrong" from "chicken is the wrong yardstick" is to redo the
#   reference-guided step against a genome that shares the grouse karyotype.
#
#   Lagopus muta (rock ptarmigan) is in Tetraoninae — the same subfamily as
#   Tympanuchus — so if these joins are real karyotype differences from chicken
#   they should become SUPPORTED, or disappear entirely, when ptarmigan is the
#   reference. If they persist against ptarmigan too, they are ours to explain.
#
# WHAT THIS DOES
#   1. Resolves the NCBI assembly directory name for the accession (it is NOT
#      hardcoded — see RESOLVING THE DIRECTORY NAME below)
#   2. Downloads genomic FASTA, genomic GFF and the assembly report
#   3. Builds the RefSeq-accession -> chromosome-name map that step 05 needs to
#      rename RagTag output to chr_1 .. chr_Z, chr_W, chr_MT
#   4. Sanity-checks the result and prints the karyotype it found, so you can
#      compare chromosome count and Z/W presence against chicken before
#      committing 46 RagTag runs to it
#
# RESOLVING THE DIRECTORY NAME
#   NCBI's FTP layout is
#     genomes/all/GCF/023/343/835/<GCF_023343835.1_ASMNAME>/
#   where ASMNAME is the submitter's assembly name. An earlier version of this
#   pipeline
#   hardcodes that directory because GRCg7b's name was known. Here the script
#   reads the parent directory listing and discovers ASMNAME itself. That is
#   deliberate: a wrong hardcoded assembly name produces a 404 at download
#   time, and the cheapest place to find out is here rather than inside an
#   array task. Set ASM_DIR_OVERRIDE below if the compute node cannot reach
#   NCBI and you want to supply the name by hand.
#
# RUNTIME: minutes, and it is a download — run it on a login node:
#     bash 03_ptarmigan_reference.sh
#
#   or as a job:
#     sbatch 03_ptarmigan_reference.sh
#
# OUTPUT, under ${PROJECT_DIR}/ref:
#   <ASM_DIR>_genomic.fna            reference FASTA (RagTag + Liftoff target)
#   <ASM_DIR>_genomic.gff            RefSeq annotation (Liftoff source)
#   <ASM_DIR>_assembly_report.txt    NCBI assembly report
#   <ASM_DIR>.chr_map.tsv            accession -> chromosome name
#   ptarmigan_asm_dir.txt            the resolved directory name, for steps 09-13
# =============================================================================
#SBATCH --job-name=ptarmigan_ref_dl
#SBATCH --output=logs/%x_%j.out
#SBATCH --error=logs/%x_%j.err
#SBATCH -A dewoody
#SBATCH -t 02:00:00
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=2
#SBATCH --mem=8G
#SBATCH -p cpu
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=blackan@purdue.edu

set -euo pipefail

module unload anaconda 2>/dev/null || true
ml biocontainers
ml samtools/1.22.1

# unset LD_PRELOAD: RCAC's XALT usage-tracking library is injected via
# LD_PRELOAD and fails on some nodes (GLIBC mismatch), which can kill subshells
# under `set -e`. It's accounting only — safe to drop.
unset LD_PRELOAD

# =============================================================================
# USER-DEFINED VARIABLES
# =============================================================================
PROJECT_DIR="${CLUSTER_SCRATCH}/GROUSE/grouse_asm"
REF_DIR="${PROJECT_DIR}/ref"

# Rock ptarmigan, assembly bLagMut1 (RefSeq). GCF_ rather than GCA_ matters
# here: RefSeq genomes are the ones that carry an NCBI annotation, which is
# what Liftoff needs as its source GFF.
ACCESSION="GCF_023343835.1"

# Leave empty to resolve the assembly directory name from NCBI (preferred).
# Set it only if the node has no outbound access and you are supplying the
# name by hand, e.g. ASM_DIR_OVERRIDE="GCF_023343835.1_bLagMut1_primary".
ASM_DIR_OVERRIDE=""

# Written so steps 09-13 can read the resolved name instead of each repeating
# the lookup.
ASM_DIR_RECORD="${REF_DIR}/ptarmigan_asm_dir.txt"

# Prefix stripped from the assembly report's Sequence-Name column when building
# the chromosome map. bLagMut1 is a Sanger/VGP curated assembly and names its
# chromosomes SUPER_1 ... SUPER_38, SUPER_Z, SUPER_W; without stripping, step 05
# would emit chr_SUPER_1 and nothing downstream would recognise it. Covers the
# other conventions that turn up too (chr1, Chr1, LG1). Set to "" to disable.
NAME_STRIP_REGEX="^(SUPER|super|Super|SCAFFOLD|CHR|Chr|chr|LG|lg)[_-]?"

# true  : abort if NCBI has no annotation GFF for this accession. Liftoff is
#         half the point of the comparison, so the default is to stop and tell
#         you rather than quietly produce an unannotated run.
# false : download what exists and let step 05 skip Liftoff.
REQUIRE_ANNOTATION=true

mkdir -p logs "$REF_DIR"

for BIN in wget samtools awk; do
    if ! command -v "$BIN" >/dev/null 2>&1; then
        echo "ERROR: '${BIN}' not on PATH."
        exit 1
    fi
done

# =============================================================================
# RESOLVE THE ASSEMBLY DIRECTORY NAME
# =============================================================================
# Builds the parent URL from the accession digits the same way step 02 does,
# then reads the listing. Every stage swallows failure explicitly: a failing
# command substitution in an assignment is fatal under `set -e` + `pipefail`,
# which is exactly how an early version of the BUSCO step in this pipeline died.
ncbi_parent_url() {
    local acc="$1" prefix digits
    prefix="${acc%%_*}"                  # GCF
    digits="${acc#*_}"; digits="${digits%%.*}"   # 023343835
    if [[ ! "$digits" =~ ^[0-9]{9}$ ]]; then
        printf ''
        return 0
    fi
    printf 'https://ftp.ncbi.nlm.nih.gov/genomes/all/%s/%s/%s/%s/' \
        "$prefix" "${digits:0:3}" "${digits:3:3}" "${digits:6:3}"
    return 0
}

resolve_asm_dir() {
    local acc="$1" url listing acc_re
    url="$(ncbi_parent_url "$acc")"
    if [[ -z "$url" ]]; then
        printf ''
        return 0
    fi
    listing="$(wget -q -O - "$url" 2>/dev/null || true)"
    if [[ -z "$listing" ]]; then
        printf ''
        return 0
    fi
    acc_re="${acc//./\\.}"               # the dot in ".1" is a regex metachar
    # Directory entries appear as href="GCF_023343835.1_NAME/". Strip anything
    # that is a file rather than the assembly directory.
    printf '%s' "$listing" \
        | grep -oE "${acc_re}_[A-Za-z0-9._-]+" \
        | grep -vE '\.(txt|gz|gbff|fna|gff|md5|xml|json)$' \
        | sort -u \
        | head -n1 || true
    return 0
}

echo ">>> Accession: ${ACCESSION}"

if [[ -n "$ASM_DIR_OVERRIDE" ]]; then
    ASM_DIR="$ASM_DIR_OVERRIDE"
    echo ">>> Assembly directory (manual override): ${ASM_DIR}"
elif [[ -s "$ASM_DIR_RECORD" ]]; then
    ASM_DIR="$(head -n1 "$ASM_DIR_RECORD")"
    echo ">>> Assembly directory (cached): ${ASM_DIR}"
else
    echo ">>> Resolving assembly directory name from NCBI"
    PARENT_URL="$(ncbi_parent_url "$ACCESSION")"
    echo "    listing: ${PARENT_URL}"
    ASM_DIR="$(resolve_asm_dir "$ACCESSION")"
    if [[ -z "$ASM_DIR" ]]; then
        echo "ERROR: could not resolve the assembly directory name for ${ACCESSION}."
        echo ""
        echo "  Either this node cannot reach ftp.ncbi.nlm.nih.gov, or the"
        echo "  accession is wrong. To check by hand, open:"
        echo "    ${PARENT_URL}"
        echo "  and look for the single ${ACCESSION}_<name> directory, then set"
        echo "    ASM_DIR_OVERRIDE=\"${ACCESSION}_<name>\""
        echo "  near the top of this script and rerun."
        exit 1
    fi
    echo ">>> Assembly directory: ${ASM_DIR}"
fi

# Guard against a resolved name that is not actually for this accession — a
# mangled listing parse would otherwise send every later step to a 404.
if [[ "$ASM_DIR" != "${ACCESSION}_"* ]]; then
    echo "ERROR: resolved directory name does not start with ${ACCESSION}_:"
    echo "       '${ASM_DIR}'"
    exit 1
fi

BASE_URL="$(ncbi_parent_url "$ACCESSION")${ASM_DIR}"

REF_FASTA="${REF_DIR}/${ASM_DIR}_genomic.fna"
REF_GFF="${REF_DIR}/${ASM_DIR}_genomic.gff"
REF_REPORT="${REF_DIR}/${ASM_DIR}_assembly_report.txt"
REF_CHR_MAP="${REF_DIR}/${ASM_DIR}.chr_map.tsv"

# =============================================================================
# DOWNLOAD
# =============================================================================
# Each file is fetched to a .part and only moved into place once wget succeeds,
# so an interrupted download can never be mistaken for a complete one by the
# `-s` guards in steps 09-13.
fetch_gz() {
    local url="$1" dest="$2" label="$3"
    if [[ -s "$dest" ]]; then
        echo "  ${label}: already present — skipping"
        return 0
    fi
    echo "  ${label}: ${url}"
    if ! wget -q -O "${dest}.gz.part" "$url"; then
        rm -f "${dest}.gz.part"
        return 1
    fi
    mv "${dest}.gz.part" "${dest}.gz"
    gunzip -f "${dest}.gz"
    if [[ ! -s "$dest" ]]; then
        echo "    ERROR: decompressed to nothing: ${dest}"
        return 1
    fi
    return 0
}

echo ""
echo ">>> Downloading into ${REF_DIR}"

if ! fetch_gz "${BASE_URL}/${ASM_DIR}_genomic.fna.gz" "$REF_FASTA" "FASTA"; then
    echo "ERROR: failed to download the genomic FASTA."
    echo "       Check ${BASE_URL}/"
    exit 1
fi

if [[ ! -s "$REF_REPORT" ]]; then
    echo "  assembly report: ${BASE_URL}/${ASM_DIR}_assembly_report.txt"
    if ! wget -q -O "${REF_REPORT}.part" "${BASE_URL}/${ASM_DIR}_assembly_report.txt"; then
        rm -f "${REF_REPORT}.part"
        echo "ERROR: failed to download the assembly report."
        exit 1
    fi
    mv "${REF_REPORT}.part" "$REF_REPORT"
else
    echo "  assembly report: already present — skipping"
fi

# The annotation is the one piece that is not guaranteed: a RefSeq genome is
# normally annotated, but not every GCF_ accession has a GFF, and finding that
# out here is far cheaper than finding out in array task 14.
ANNOTATION_PRESENT=true
if ! fetch_gz "${BASE_URL}/${ASM_DIR}_genomic.gff.gz" "$REF_GFF" "GFF annotation"; then
    ANNOTATION_PRESENT=false
    rm -f "$REF_GFF"
    echo ""
    echo "  !! No annotation GFF at ${BASE_URL}/${ASM_DIR}_genomic.gff.gz"
    echo "     This accession appears not to carry an NCBI annotation."
    if [[ "$REQUIRE_ANNOTATION" == true ]]; then
        echo ""
        echo "ERROR: REQUIRE_ANNOTATION=true and no GFF is available."
        echo "       Options:"
        echo "         - check ${BASE_URL}/ for a differently named GFF"
        echo "         - pick an annotated Tetraoninae/Phasianidae assembly instead"
        echo "         - set REQUIRE_ANNOTATION=false to run RagTag + BUSCO only"
        echo "           (step 05 will then skip Liftoff, and the karyotype and"
        echo "           join-support comparisons still work — only annotation"
        echo "           transfer is lost)"
        exit 1
    fi
    echo "     REQUIRE_ANNOTATION=false — continuing without it; step 05 will"
    echo "     skip Liftoff."
fi

# =============================================================================
# CHROMOSOME-NAME MAP
# Identical construction to steps 02 and 05 so chr_* names mean the same
# thing in both trees: RefSeq accession -> chromosome name, assembled molecules
# only. Report columns: 1 Sequence-Name, 2 Sequence-Role, ..., 7 RefSeq-Accn.
# =============================================================================
echo ""
echo ">>> Building chromosome-name map"
if [[ ! -s "$REF_CHR_MAP" ]]; then
    # Sequence-Name is whatever the submitter used. Sanger/VGP curated assemblies
    # (which this is) name curated scaffolds SUPER_1, SUPER_2, ..., SUPER_Z,
    # SUPER_W. Left alone, step 05 would prefix "chr_" and produce chr_SUPER_1,
    # which (a) is invisible to step 07's EXTRACT_SEQS=(chr_W chr_Z chr_MT), so
    # the sex and organelle sequences would never be split out, and (b) fails the
    # chr_(\d+) regex in 09_busco_plots.py's natural_chr_key, so every chromosome
    # sorts into the unplaced-scaffold bucket and the karyotype ordering silently
    # collapses. Strip the prefix here, once, so chr_1 / chr_Z / chr_W / chr_MT
    # are what every downstream step sees.
    awk -F'\t' -v re="$NAME_STRIP_REGEX" '
        !/^#/ && $2 == "assembled-molecule" {
            name = $1
            if (re != "") sub(re, "", name)
            print $7"\t"name
        }' "$REF_REPORT" > "${REF_CHR_MAP}.part"
    if [[ ! -s "${REF_CHR_MAP}.part" ]]; then
        echo "ERROR: chromosome map came out empty — check ${REF_REPORT}"
        rm -f "${REF_CHR_MAP}.part"
        exit 1
    fi
    mv "${REF_CHR_MAP}.part" "$REF_CHR_MAP"

    # Show what the stripping did, so a surprising convention is visible rather
    # than applied silently.
    N_STRIPPED=$(awk -F'\t' -v re="$NAME_STRIP_REGEX" \
        '!/^#/ && $2 == "assembled-molecule" && $1 ~ re {n++} END {print n+0}' "$REF_REPORT")
    if (( N_STRIPPED > 0 )); then
        echo "  stripped a name prefix from ${N_STRIPPED} sequence name(s), e.g.:"
        awk -F'\t' -v re="$NAME_STRIP_REGEX" \
            '!/^#/ && $2 == "assembled-molecule" && $1 ~ re {
                 name=$1; sub(re,"",name); printf "       %-14s -> chr_%s\n", $1, name
             }' "$REF_REPORT" | head -4
    fi
fi
N_CHR=$(wc -l < "$REF_CHR_MAP")
echo "  ${REF_CHR_MAP}: ${N_CHR} assembled molecules"

# Steps 05 and 08 build names as "chr_" + column 2. If the submitter already named
# sequences "chr1" or "1chr", that would yield chr_chr1 and break every
# downstream name comparison against the chicken tree. Catch it here.
if awk -F'\t' '$2 ~ /^[Cc][Hh][Rr]/ {found=1} END {exit !found}' "$REF_CHR_MAP"; then
    echo ""
    echo "  !! WARNING: some sequence names STILL begin with 'chr' after stripping:"
    awk -F'\t' '$2 ~ /^[Cc][Hh][Rr]/ {printf "       %s -> %s\n", $1, $2}' "$REF_CHR_MAP" | head -5
    echo "     Step 05 prefixes 'chr_', which would give names like chr_chr1."
    echo "     Extend NAME_STRIP_REGEX at the top of this script, or edit column 2"
    echo "     of ${REF_CHR_MAP}, before running step 05."
fi

# =============================================================================
# INDEX AND REPORT THE KARYOTYPE
# =============================================================================
if [[ ! -s "${REF_FASTA}.fai" || "${REF_FASTA}.fai" -ot "$REF_FASTA" ]]; then
    echo ""
    echo ">>> Indexing FASTA"
    samtools faidx "$REF_FASTA"
fi

echo ""
echo "============================================================"
echo ">>> Reference karyotype as NCBI reports it"
echo "============================================================"
# Join the map onto the .fai so lengths are shown against chromosome names.
# Sorted by length descending, which is how RagTag output gets numbered.
awk -F'\t' -v map="$REF_CHR_MAP" '
    BEGIN {
        while ((getline line < map) > 0) {
            split(line, a, "\t")
            chrname[a[1]] = a[2]
        }
    }
    {
        if ($1 in chrname) printf "%s\tchr_%s\t%d\n", $1, chrname[$1], $2
    }
' "${REF_FASTA}.fai" \
    | sort -k3,3nr \
    | awk -F'\t' '{printf "  %-16s %-10s %12.2f Mb\n", $1, $2, $3/1e6}'

# =============================================================================
# IS THE NUMBERING A SIZE RANK, OR A HOMOLOGY ASSIGNMENT?
# This decides what every downstream chromosome name MEANS, so it is checked
# rather than assumed.
#
# Two conventions are in use. Some references number chromosomes by homology to
# an established karyotype, so chrN is "the chromosome everyone calls N". Curated
# assemblies out of the Sanger/VGP pipeline instead number scaffolds by DESCENDING
# LENGTH, so chrN is merely "the Nth longest" and carries no homology claim at all.
#
# Confusing the two is silent and serious: scaffolding the grouse against a
# size-ranked reference and then calling the result chr_6 produces a "chr_6" that
# is a different chromosome from chr_6 in the chicken-based tree, and every
# cross-tree comparison inherits the error without any step failing.
#
# The tell is monotonicity. Under homology naming, sizes wander (chicken chr4 is
# 91 Mb, chr5 59 Mb, chr6 36 Mb — descending, but chicken is the karyotype others
# are named after). Under size-rank naming, sizes descend essentially perfectly
# AND the numbering steps over the lettered sex chromosomes.
# =============================================================================
echo ""
echo "============================================================"
echo ">>> Does the numbering encode homology, or just size rank?"
echo "============================================================"
awk -F'\t' -v map="$REF_CHR_MAP" '
    BEGIN {
        while ((getline line < map) > 0) { split(line, a, "\t"); nm[a[1]] = a[2] }
    }
    { if (($1 in nm) && nm[$1] ~ /^[0-9]+$/) len[nm[$1] + 0] = $2 }
    END {
        n = 0
        for (i = 1; i <= 200; i++) if (i in len) { n++; idx[n] = i }
        if (n < 5) { print "  too few numbered chromosomes to judge"; exit }
        viol = 0
        for (j = 2; j <= n; j++) if (len[idx[j]] > len[idx[j-1]]) viol++
        pairs = n - 1
        printf "  numbered chromosomes            : %d\n", n
        printf "  consecutive pairs in size order : %d / %d\n", pairs - viol, pairs
        if (viol <= pairs * 0.1) {
            print ""
            print "  !! SIZE-RANK NAMING DETECTED."
            print "     The numbers are a length ordering, NOT homology to chicken."
            print "     Do NOT read this reference chrN as chicken chrN, and do not"
            print "     let step 05 name grouse chromosomes from it until the"
            print "     correspondence has been computed:"
            print ""
            print "         sbatch 04_chromosome_homology.sh"
            print ""
            print "     That aligns this reference to chicken and writes a homology"
            print "     map, which step 05 can then use (CHR_NAMING=homology) so chr_N"
            print "     means the same chromosome in both trees."
        } else {
            printf "  -> %d of %d pairs out of size order: numbering looks homology-based.\n", viol, pairs
            print "     Running 04_chromosome_homology.sh is still worthwhile to confirm."
        }
    }
' "${REF_FASTA}.fai"

TOTAL_BP=$(awk -F'\t' '{s+=$2} END {print s+0}' "${REF_FASTA}.fai")
N_SEQ=$(wc -l < "${REF_FASTA}.fai")
HAS_Z=$(awk -F'\t' '$2 == "Z" {print "yes"; exit}' "$REF_CHR_MAP")
HAS_W=$(awk -F'\t' '$2 == "W" {print "yes"; exit}' "$REF_CHR_MAP")
HAS_MT=$(awk -F'\t' '$2 == "MT" {print "yes"; exit}' "$REF_CHR_MAP")

echo ""
echo "  total sequences    : ${N_SEQ}"
printf "  total length       : %.2f Gb\n" "$(awk -v b="$TOTAL_BP" 'BEGIN{print b/1e9}')"
echo "  assembled molecules: ${N_CHR}"
echo "  chr_Z present      : ${HAS_Z:-no}"
echo "  chr_W present      : ${HAS_W:-no}"
echo "  chr_MT present     : ${HAS_MT:-no}"
echo ""
echo "  For contrast, chicken GRCg7b has 42 assembled molecules"
echo "  (chr_1..chr_39 plus chr_W, chr_Z, chr_MT) and ~1.05 Gb."
echo ""
echo "  WORTH A LOOK BEFORE COMMITTING 46 RagTag RUNS:"
echo "    - Does the autosome count match chicken's 39? A difference is a"
echo "      karyotype difference, a missing microchromosome, or both."
echo "    - Is there a chromosome near 65.8 Mb? That is the size of chicken"
echo "      chr6 + chr8 combined, and of the single chromosome seen in all 46"
echo "      grouse haplotypes. A match is a CANDIDATE for the shared fusion."
echo "    - Is there a chromosome matching chicken chr2 (~149 Mb)? If not, that"
echo "      chromosome is split here, which is a second karyotype difference."
echo ""
echo "    Size matching is a hypothesis generator ONLY. Two chromosomes can"
echo "    share a length without sharing ancestry. 04_chromosome_homology.sh"
echo "    settles it by alignment; do not draw conclusions from the sizes alone."

# Record the resolved name for steps 09-13.
printf '%s\n' "$ASM_DIR" > "$ASM_DIR_RECORD"

echo ""
echo "============================================================"
echo ">>> Complete: $(date)"
echo "  FASTA       : ${REF_FASTA}"
if [[ "$ANNOTATION_PRESENT" == true ]]; then
    echo "  GFF         : ${REF_GFF}"
else
    echo "  GFF         : NOT AVAILABLE — step 05 will skip Liftoff"
fi
echo "  report      : ${REF_REPORT}"
echo "  chr map     : ${REF_CHR_MAP}"
echo "  dir record  : ${ASM_DIR_RECORD}"
echo ""
echo "  Next:"
echo "    N=\$(grep -v '^#' assembly_manifest.tsv | tail -n +2 | grep -c .)"
echo "    sbatch --array=0-\$((N-1))%6 05_ragtag_liftoff_array.sh"
echo "============================================================"
