#!/bin/bash
# =============================================================================
# DOWNLOAD THE CHICKEN REFERENCE (GCF_016699485.2 / GRCg7b)
# Step 01 — run once. Needed only by 04_chromosome_homology.sh.
#
# WHY THE CHICKEN GENOME IS STILL HERE
#   These assemblies are ordered against the rock ptarmigan (step 03), not
#   against chicken. Chicken is downloaded for exactly one purpose: to work out
#   what each ptarmigan chromosome number MEANS.
#
#   bLagMut1 is a Sanger/VGP curated assembly, and it numbers its chromosomes
#   SUPER_1, SUPER_2, ... in DESCENDING LENGTH ORDER. That is a size rank, not a
#   homology statement: its SUPER_6 is not chicken chromosome 6. Every published
#   galliform karyotype, gene name and cytogenetic result, on the other hand, is
#   expressed in chicken (GGA) coordinates. So a ptarmigan chromosome number is
#   unusable on its own — it has to be tied to the chicken karyotype by
#   alignment before it can be related to anything in the literature.
#
#   Step 04 does that alignment. This script supplies its target genome and the
#   accession -> chromosome-name map it joins on. Nothing else in the pipeline
#   reads the chicken genome: RagTag and Liftoff both use ptarmigan.
#
#   The FASTA only. No GFF is downloaded, because no annotation is transferred
#   from chicken any more — Liftoff (step 05) uses the ptarmigan RefSeq
#   annotation, which is the whole point of having changed reference.
#
# OUTPUT, under ${PROJECT_DIR}/ref:
#   GCF_016699485.2_bGalGal1.mat.broiler.GRCg7b_genomic.fna      + .fai
#   GCF_016699485.2_bGalGal1.mat.broiler.GRCg7b_assembly_report.txt
#   GCF_016699485.2_bGalGal1.mat.broiler.GRCg7b.chr_map.tsv
#
# Idempotent: anything already present is left alone.
#
#   sbatch 01_download_chicken_reference.sh
#
# or, since it is a download and two awk passes, just:
#
#   bash 01_download_chicken_reference.sh
# =============================================================================
#SBATCH --job-name=chicken_ref
#SBATCH --output=logs/%x_%j.out
#SBATCH --error=logs/%x_%j.err
#SBATCH -A fnrdewoody
#SBATCH -t 0-04:00:00
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=8G
#SBATCH -p cpu
#SBATCH --mail-type=BEGIN,END,FAIL
#SBATCH --mail-user=blackan@purdue.edu

set -euo pipefail

# unset LD_PRELOAD: RCAC's XALT library is injected this way and fails on some
# nodes (GLIBC mismatch), which can kill subshells under `set -e`.
unset LD_PRELOAD

# =============================================================================
# USER-DEFINED VARIABLES
# =============================================================================
PROJECT_DIR="${CLUSTER_SCRATCH}/GROUSE/grouse_asm"
REF_DIR="${PROJECT_DIR}/ref"

# Hardcoded, unlike the ptarmigan accession in step 03. GRCg7b is a long-stable
# assembly and the directory name is not going to move under us; step 03 resolves
# its directory at runtime because that accession was new at the time.
ASM_DIR="GCF_016699485.2_bGalGal1.mat.broiler.GRCg7b"
BASE_URL="https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/016/699/485/${ASM_DIR}"

REF_FASTA="${REF_DIR}/${ASM_DIR}_genomic.fna"
REF_REPORT="${REF_DIR}/${ASM_DIR}_assembly_report.txt"
REF_CHR_MAP="${REF_DIR}/${ASM_DIR}.chr_map.tsv"

# Same safety net as step 03. GRCg7b names its assembled molecules 1..39, Z, W,
# MT, so nothing should match — but if NCBI ever re-releases it with a prefix,
# silently producing chr_CHR1 downstream would be worse than stripping here.
NAME_STRIP_REGEX="^(SUPER|super|Super|SCAFFOLD|CHR|Chr|chr|LG|lg)[_-]?"

mkdir -p logs "$REF_DIR"

echo ">>> 01_download_chicken_reference.sh"
echo ">>> Reference : ${ASM_DIR}"
echo ">>> Target dir: ${REF_DIR}"
echo ">>> Start     : $(date)"

for BIN in wget gunzip samtools awk; do
    if ! command -v "$BIN" >/dev/null 2>&1; then
        echo "ERROR: '${BIN}' not on PATH."
        exit 1
    fi
done

# =============================================================================
# FASTA
# Downloaded to .part and only then moved into place, so an interrupted
# download never leaves a truncated FASTA that later steps happily index.
# =============================================================================
echo ""
if [[ -s "$REF_FASTA" ]]; then
    echo ">>> FASTA already present — skipping download"
else
    echo ">>> Downloading FASTA: ${BASE_URL}/${ASM_DIR}_genomic.fna.gz"
    if ! wget -q -O "${REF_FASTA}.gz.part" "${BASE_URL}/${ASM_DIR}_genomic.fna.gz"; then
        rm -f "${REF_FASTA}.gz.part"
        echo "ERROR: download failed. Check ${BASE_URL}/"
        exit 1
    fi
    mv "${REF_FASTA}.gz.part" "${REF_FASTA}.gz"
    gunzip -f "${REF_FASTA}.gz"
    echo "    ${REF_FASTA}"
fi

# =============================================================================
# ASSEMBLY REPORT
# =============================================================================
echo ""
if [[ -s "$REF_REPORT" ]]; then
    echo ">>> Assembly report already present — skipping download"
else
    echo ">>> Downloading assembly report"
    if ! wget -q -O "${REF_REPORT}.part" "${BASE_URL}/${ASM_DIR}_assembly_report.txt"; then
        rm -f "${REF_REPORT}.part"
        echo "ERROR: failed to download the assembly report."
        exit 1
    fi
    mv "${REF_REPORT}.part" "$REF_REPORT"
    echo "    ${REF_REPORT}"
fi

# =============================================================================
# CHROMOSOME-NAME MAP
# Identical construction to step 03, so the two maps can be joined on without
# any per-reference special-casing in step 04. Report columns:
#   1 Sequence-Name, 2 Sequence-Role, ..., 7 RefSeq-Accn
# Assembled molecules only: unplaced and unlocalised scaffolds are deliberately
# left out, because a chromosome correspondence computed against them is noise.
# =============================================================================
echo ""
echo ">>> Building chromosome-name map"
if [[ -s "$REF_CHR_MAP" ]]; then
    echo "    already present — skipping"
else
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
fi
echo "    ${REF_CHR_MAP}: $(wc -l < "$REF_CHR_MAP") assembled molecules"

# Step 04 and step 05 both build names as "chr_" + column 2, so a name that
# already begins with "chr" would produce chr_chr1 and break every name
# comparison. Catch it here rather than three steps downstream.
if awk -F'\t' '$2 ~ /^[Cc][Hh][Rr]/ {found=1} END {exit !found}' "$REF_CHR_MAP"; then
    echo ""
    echo "  !! WARNING: some sequence names STILL begin with 'chr' after stripping:"
    awk -F'\t' '$2 ~ /^[Cc][Hh][Rr]/ {printf "       %s -> %s\n", $1, $2}' \
        "$REF_CHR_MAP" | head -5
    echo "     Extend NAME_STRIP_REGEX at the top of this script, or edit column 2"
    echo "     of ${REF_CHR_MAP}, before running step 04."
fi

# =============================================================================
# INDEX
# =============================================================================
if [[ ! -s "${REF_FASTA}.fai" || "${REF_FASTA}.fai" -ot "$REF_FASTA" ]]; then
    echo ""
    echo ">>> Indexing FASTA"
    samtools faidx "$REF_FASTA"
fi

# =============================================================================
# REPORT THE KARYOTYPE
# One line per assembled molecule, longest first. GRCg7b should show 1..39 plus
# Z, W and MT; chromosome 1 near 197 Mb and the microchromosomes trailing off
# below 1 Mb. If this looks nothing like that, the download is wrong.
# =============================================================================
echo ""
echo "============================================================"
echo ">>> Chicken karyotype as downloaded"
echo "============================================================"
awk -F'\t' 'NR==FNR {name[$1]=$2; next}
            ($1 in name) {printf "  chr_%-5s %12.2f Mb\n", name[$1], $2/1e6}' \
    "$REF_CHR_MAP" "${REF_FASTA}.fai" \
    | sort -k2,2 -nr | head -12
echo "  ..."
echo ""
echo ">>> Step 01 complete — $(date)"
echo ""
echo "  Next:"
echo "    bash   03_ptarmigan_reference.sh    # the working reference"
echo "    sbatch 04_chromosome_homology.sh    # needs both references"
echo "============================================================"
