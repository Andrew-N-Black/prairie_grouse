#!/bin/bash
# =============================================================================
# NCBI SUBMISSION PREP: 46 haplotype assemblies (23 birds) + mitogenomes
# Builds the files for NCBI's diploid genome submission portal, following NCBI's
# instructions for this project: pairs per bird; the Z, W and mitogenome in the
# haplotype submitted as "Haplotype 1"; no RagTag reference accessions in any
# sequence identifier; chromosome and organelle assignments in definition lines.
#
# STAGES (run from the folder holding this script and ncbi_submission_prep.py)
#   bash ncbi_submission_prep.sh manifest
#       login node, ~10 min. Writes ${OUT}/birds.tsv (species, sample, sex from W
#       content, both haplotype FASTAs, the mitogenome found under MITO_DIR) and
#       prints the sequence-name patterns of one assembly. CHECK birds.tsv: every
#       bird needs the right mitogenome FASTA (one sequence) and the right sex.
#   ONLY=F5503 sbatch --export=ALL,STAGE=prep ncbi_submission_prep.sh   # one-bird test
#   sbatch --export=ALL,STAGE=prep ncbi_submission_prep.sh                # all 23 birds
#       Writes ${OUT}/files/<sample>_hap{1,2}.fsa.gz (submission labels), rename
#       maps, haplotype_table.tsv (one row per assembly, for the portal's genome
#       info sheet) and prep_report.txt (every moved, dropped or unplaced piece).
#   sbatch --export=ALL,STAGE=fcs ncbi_submission_prep.sh                 # optional
#       NCBI FCS-adaptor and FCS-GX on every file, if FCS_DIR and GXDB are set.
#       Reports only: review ${OUT}/fcs/*/ and clean with fcs.py if anything is flagged.
#
# Female sex chromosomes: the haplotype with more chr_Z becomes Haplotype 1. Z/W
# pieces in the other haplotype are moved into Haplotype 1 only if < DUP_COV (default 20%) of the
# piece aligns to Haplotype 1 Z/W; otherwise it is treated as (partly) duplicated and
# left unplaced in Haplotype 2). Needs minimap2; without it nothing is moved.
# =============================================================================
#SBATCH --job-name=ncbi_prep
#SBATCH -A dewoody
#SBATCH -p cpu
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=16
#SBATCH --mem=48G
#SBATCH -t 12:00:00
#SBATCH --output=logs/%x_%j.out
#SBATCH --error=logs/%x_%j.err
set -euo pipefail
unset LD_PRELOAD

PROJECT_DIR="${CLUSTER_SCRATCH}/GROUSE/grouse_asm"
ASSEMBLY_DIR="${ASSEMBLY_DIR:-${PROJECT_DIR}/final_ptarmigan}"     # <SP>_<sample>_hap{1,2}.pseudo_chr.fasta
SUFFIX="${SUFFIX:-.pseudo_chr.fasta}"
MITO_DIR="${MITO_DIR:-${PROJECT_DIR}/mitogenomes}"                 # EDIT: where the 23 final mitogenomes are
OUT="${OUT:-${PROJECT_DIR}/ncbi_submission}"
ONLY="${ONLY:-}"
THREADS="${SLURM_CPUS_PER_TASK:-8}"
STAGE="${STAGE:-${1:-}}"
HERE="${SLURM_SUBMIT_DIR:-$(pwd)}"
# FCS (optional): folder with run_fcsadaptor.sh, fcs.py, fcs-adaptor.sif, fcs-gx.sif; GX database path
FCS_DIR="${FCS_DIR:-}"
GXDB="${GXDB:-}"

mkdir -p logs "$OUT"

case "$STAGE" in
manifest)
    python3 "${HERE}/ncbi_submission_prep.py" manifest --assembly-dir "$ASSEMBLY_DIR" --suffix "$SUFFIX" \
        --mito-dir "$MITO_DIR" --out "$OUT"
    ;;
prep)
    [[ -s "${OUT}/birds.tsv" ]] || { echo "ERROR: run the manifest stage first"; exit 1; }
    module load biocontainers 2>/dev/null || true
    module load minimap2 2>/dev/null || true
    command -v minimap2 >/dev/null || echo "WARNING: minimap2 not found; female Z/W pieces will not be moved"
    python3 "${HERE}/ncbi_submission_prep.py" prep --out "$OUT" --threads "$THREADS" --only "$ONLY" --dup-cov "${DUP_COV:-0.2}"
    ( cd "${OUT}/files" && md5sum *.fsa.gz > md5sums.txt )
    cat > "${OUT}/submission_note.txt" << 'EOF'
Note for the NCBI diploid genome submission (paste into the comments):

These are 46 haplotype-resolved assemblies from 23 birds of three Tympanuchus species, arranged in pairs by bird.
As advised by NCBI (genomes@ncbi.nlm.nih.gov), please place the 23 per-bird umbrella BioProjects under one
overarching umbrella BioProject for the whole study [we request that this be PRJNA1513026 / a new umbrella].
We have chosen to include the Z and W chromosomes and the mitochondrial genome in Haplotype 1 of each bird.
In the six females, Haplotype 1 is the haplotype carrying the Z and W; the second haplotype has no Z or W.
Sequence identifiers contain no accessions from the reference (rock ptarmigan, GCF_023343835.1) used for
RagTag ordering. Chromosome and organelle assignments are given in the FASTA definition lines.
EOF
    echo "Files: ${OUT}/files   Table: ${OUT}/haplotype_table.tsv   Report: ${OUT}/prep_report.txt"
    ;;
fcs)
    [[ -n "$FCS_DIR" && -n "$GXDB" ]] || { echo "ERROR: set FCS_DIR and GXDB"; exit 1; }
    taxid() {   # NCBI taxonomy ID for a species name (E-utilities)
        curl -s "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esearch.fcgi?db=taxonomy&term=$(echo "$1" | sed 's/ /%20/g')%5BScientific%20Name%5D" \
            | grep -o '<Id>[0-9]*</Id>' | head -1 | tr -dc '0-9'
    }
    declare -A TAX=( [LEPC]="${TAXID_LEPC:-$(taxid 'Tympanuchus pallidicinctus')}"
                     [GRPC]="${TAXID_GRPC:-$(taxid 'Tympanuchus cupido')}"
                     [STGR]="${TAXID_STGR:-$(taxid 'Tympanuchus phasianellus')}" )
    for k in LEPC GRPC STGR; do
        [[ -n "${TAX[$k]}" ]] || { echo "ERROR: no taxid for ${k}; set TAXID_${k}"; exit 1; }
        echo "taxid ${k}: ${TAX[$k]}"
    done
    export FCS_DEFAULT_IMAGE="${FCS_DIR}/fcs-gx.sif"
    tail -n +2 "${OUT}/haplotype_table.tsv" | while IFS=$'\t' read -r _ sample sp _ _ _ _ _ fn _; do
        d="${OUT}/fcs/${fn%.fsa.gz}"; mkdir -p "$d/adaptor" "$d/gx"
        echo ">>> ${fn}"
        bash "${FCS_DIR}/run_fcsadaptor.sh" --fasta-input "${OUT}/files/${fn}" --output-dir "$d/adaptor" \
            --euk --container-engine singularity --image "${FCS_DIR}/fcs-adaptor.sif"
        python3 "${FCS_DIR}/fcs.py" screen genome --fasta "${OUT}/files/${fn}" --out-dir "$d/gx" \
            --gx-db "$GXDB" --tax-id "${TAX[$sp]}"
    done
    echo "Reports in ${OUT}/fcs/<assembly>/{adaptor,gx}; anything flagged needs fcs.py clean before upload"
    ;;
*)
    echo "Usage: bash $0 manifest | sbatch --export=ALL,STAGE=prep $0 | sbatch --export=ALL,STAGE=fcs $0" >&2
    exit 1
    ;;
esac
