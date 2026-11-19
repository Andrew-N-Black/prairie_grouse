#!/bin/bash
# =============================================================================
# SLURM JOB SUBMISSION: STRUCTURAL VARIANTS THAT DIFFERENTIATE SPECIES
# Step 13 — requires step 12 (per-chromosome graphs in pangenome/<CHROM>/ and
# comparison/genome_length_by_haplotype.tsv).
#
# Answers two questions for the manuscript:
#   1. How many structural variants (SVs, >= 50 bp) are there, and what is
#      their size distribution, overall and among those that separate species?
#   2. Are species-diagnostic SVs concentrated on the Z relative to autosomes?
#
# Per chromosome (array task, STAGE=chrom):
#   a. choose a reference haplotype: the longest haplotype that carries the
#      chromosome as ONE contig (no _p pieces), taking species in REF_ORDER
#      (default LEPC first, so autosomes and the Z use a male reference; chr_W
#      falls through to a female STGR or GRPC haplotype automatically)
#   b. vg deconstruct the chromosome graph against that path (-a -e: every
#      snarl, path traversals only); genotypes are PanSN haplotypes, so each
#      bird is one VCF sample with hap1|hap2
#   c. vcfbub -l 0 -a MAX_BUBBLE: keep top-level bubbles and pop any larger
#      than MAX_BUBBLE into their children, then vcfwave -I 1000 to realign
#      each bubble into SNPs, indels and inversions (as in pggb -V)
#   d. count, per species, the haplotypes that carry the chromosome (total path
#      length >= half the chromosome's median), so a female haplotype holding
#      only a fragment of the Z is not counted as a missing Z call
# Summary (STAGE=summary): 13_sv_summary.py classifies every SV allele
#   (INS / DEL / COMPLEX) and calls it
#     - diagnostic for species X: fixed (or absent) in X and the opposite state
#       fixed in both other species
#     - a fixed difference for a species pair: fixed in one, absent in the other
#   both with all called haplotypes and after subsampling every species to
#   EQUAL_N haplotypes per site (fixation is easier to reach in a small sample,
#   and GRPC has only 7 Z haplotypes against 12 autosomal ones), then compares
#   the Z with the autosomes (Fisher's exact test on the fraction of SVs that
#   are diagnostic, and density per Mb of reference sequence).
#
# USAGE (from the repo root, on a login node)
#   bash 13_structural_variants.sh submit           # array over chromosomes + summary job
#   CHROM=chr_20 sbatch --export=ALL,STAGE=chrom 13_structural_variants.sh   # one-chromosome test
#   STAGE=summary sbatch --export=ALL -c 4 --mem=32G -t 04:00:00 13_structural_variants.sh
#
# OUTPUT: ${SV_DIR}/<CHROM>/{ref.txt,haps_per_species.tsv,<CHROM>.sv.vcf.gz}
#         ${SV_DIR}/summary/{sv_table.tsv.gz, sv_per_chromosome.tsv,
#            sv_size_distribution.tsv, sv_Z_vs_autosomes.tsv, sv_summary.txt}
#
# SOFTWARE: vg from ${CONTAINER_DIR}/vg.sif and vcfbub from pggb.sif (both pulled
#   by accessory/pull_containers.sh); the summary needs only Python 3 (stdlib).
# =============================================================================
#SBATCH --job-name=sv_species
#SBATCH -A dewoody
#SBATCH -p cpu
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=96
#SBATCH --mem=190G
#SBATCH -t 1-00:00:00
#SBATCH --output=logs/%x_%A_%a.out
#SBATCH --error=logs/%x_%A_%a.err
#SBATCH --mail-type=FAIL
#SBATCH --mail-user=blackan@purdue.edu
set -euo pipefail
unset LD_PRELOAD   # RCAC XALT preload breaks some subshells (GLIBC mismatch)

# ----------------------------- settings ---------------------------------------
PROJECT_DIR="${CLUSTER_SCRATCH}/GROUSE/grouse_asm"
PANGENOME_DIR="${PROJECT_DIR}/pangenome"
CHROM_LIST="${CHROM_LIST:-${PANGENOME_DIR}/chromosomes.txt}"
LENGTHS="${PROJECT_DIR}/comparison/genome_length_by_haplotype.tsv"
SV_DIR="${SV_DIR:-${PROJECT_DIR}/structural_variants}"
CONTAINER_DIR="${PROJECT_DIR}/containers"
VG_SIF="${VG_SIF:-${CONTAINER_DIR}/vg.sif}"
PGGB_SIF="${PGGB_SIF:-${CONTAINER_DIR}/pggb.sif}"
REF_ORDER="${REF_ORDER:-LEPC STGR GRPC}"   # species searched, in order, for the reference haplotype
MAX_BUBBLE="${MAX_BUBBLE:-100000}"          # vcfbub -a: pop bubbles larger than this (bp)
MIN_SV="${MIN_SV:-50}"                      # minimum SV length (bp)
MIN_CALLED_FRAC="${MIN_CALLED_FRAC:-0.7}"   # per species: called haplotypes / haplotypes carrying the chromosome
EQUAL_N="${EQUAL_N:-6}"                     # haplotypes per species after subsampling
SEED="${SEED:-1}"
THREADS="${SLURM_CPUS_PER_TASK:-8}"
STAGE="${STAGE:-${1:-}}"
REPO="${SLURM_SUBMIT_DIR:-$(pwd)}"          # sbatch runs a spooled copy, so locate the .py by submit dir

APPTAINER_BIN="$(command -v apptainer || command -v singularity || true)"
cx() { local sif="$1"; shift; "$APPTAINER_BIN" exec --cleanenv --bind "${CLUSTER_SCRATCH}" "$sif" "$@"; }

# ----------------------------- submit -----------------------------------------
if [[ "$STAGE" == "submit" ]]; then
    [[ -s "$CHROM_LIST" ]] || { echo "ERROR: ${CHROM_LIST} not found"; exit 1; }
    [[ -s "$LENGTHS" ]]    || { echo "ERROR: ${LENGTHS} not found (run 03 STAGE=summary)"; exit 1; }
    for f in "$VG_SIF" "$PGGB_SIF"; do [[ -s "$f" ]] || { echo "ERROR: ${f} missing (accessory/pull_containers.sh)"; exit 1; }; done
    [[ -s "${REPO}/13_sv_summary.py" ]] || { echo "ERROR: run from the repo root (13_sv_summary.py not found)"; exit 1; }
    mkdir -p logs "$SV_DIR"
    N=$(grep -c . "$CHROM_LIST")
    jid=$(sbatch --parsable --export=ALL,STAGE=chrom --array=0-$((N - 1))%10 "$0")
    sid=$(sbatch --parsable --export=ALL,STAGE=summary --dependency=afterok:"$jid" \
          --job-name=sv_summary -c 4 --mem=32G -t 04:00:00 "$0")
    echo "Submitted array ${jid} (${N} chromosomes) and summary ${sid} (runs after all tasks succeed)"
    exit 0
fi

# ----------------------------- one chromosome ---------------------------------
if [[ "$STAGE" == "chrom" ]]; then
    if [[ -z "${CHROM:-}" ]]; then
        CHROM=$(sed -n "$((SLURM_ARRAY_TASK_ID + 1))p" "$CHROM_LIST")
    fi
    [[ -n "$CHROM" ]] || { echo "ERROR: no chromosome for this task"; exit 1; }
    OUT="${SV_DIR}/${CHROM}"; mkdir -p "$OUT"
    GFA=$(find "${PANGENOME_DIR}/${CHROM}" -maxdepth 1 -name "*.final.gfa" | head -n1)
    [[ -s "$GFA" ]] || { echo "ERROR: no *.final.gfa in ${PANGENOME_DIR}/${CHROM}"; exit 1; }
    echo ">>> ${CHROM}: $(date)  graph ${GFA}"

    # a. reference haplotype + d. haplotypes per species carrying this chromosome
    python3 - "$LENGTHS" "$CHROM" "$REF_ORDER" "$OUT" << 'PYEOF'
import sys, csv, re, statistics
lengths, chrom, order, out = sys.argv[1], sys.argv[2], sys.argv[3].split(), sys.argv[4]
tot, pieces = {}, {}
for r in csv.DictReader(open(lengths), delimiter="\t"):
    c = r["contig"]
    if re.sub(r"_p\d+$", "", c) != chrom:
        continue
    k = (r["species"], r["sample"], r["haplotype"])
    tot[k] = tot.get(k, 0) + int(r["path_length_bp"])
    pieces.setdefault(k, []).append(c)
if not tot:
    sys.exit(f"ERROR: {chrom} not in {lengths}")
med = statistics.median(tot.values())
carry = {k: v for k, v in tot.items() if v >= 0.5 * med}
with open(f"{out}/haps_per_species.tsv", "w") as f:
    f.write("species\tn_haplotypes\tchrom_median_bp\n")
    for sp in sorted({k[0] for k in tot}):
        f.write(f"{sp}\t{sum(1 for k in carry if k[0] == sp)}\t{int(med)}\n")
ref = None
for sp in order:
    cand = [(v, k) for k, v in carry.items() if k[0] == sp and pieces[k] == [chrom]]
    if cand:
        v, k = max(cand)
        ref = (f"{k[0]}.{k[1]}#{k[2].replace('hap', '')}#{chrom}", v)
        break
if ref is None:
    sys.exit(f"ERROR: no single-contig haplotype for {chrom}")
open(f"{out}/ref.txt", "w").write(f"{ref[0]}\t{ref[1]}\n")
print(f"  reference path {ref[0]} ({ref[1]:,} bp); haplotypes carrying {chrom}: "
      + ", ".join(f"{sp} {sum(1 for k in carry if k[0] == sp)}" for sp in sorted({k[0] for k in tot})))
PYEOF
    REF=$(cut -f1 "${OUT}/ref.txt")

    # b. deconstruct (older vg needs -H '#' to split PanSN names; newer vg parses them itself)
    HFLAG=()
    cx "$VG_SIF" vg deconstruct --help 2>&1 | grep -q -- '--path-sep\|-H,' && HFLAG=(-H '#')
    echo "  vg deconstruct -P ${REF} ${HFLAG[*]:-} -a -e"
    cx "$VG_SIF" vg deconstruct -P "$REF" "${HFLAG[@]}" -a -e -t "$THREADS" "$GFA" > "${OUT}/${CHROM}.raw.vcf"

    # c. pop large bubbles, then e. decompose each bubble into its primitive variants with vcfwave
    #    (as pggb -V does). Without this, a long bubble that differs only by a few SNPs looks like a
    #    "complex" SV, and a bubble with many haplotype-specific alleles can never look fixed.
    cx "$PGGB_SIF" vcfbub -l 0 -a "$MAX_BUBBLE" --input "${OUT}/${CHROM}.raw.vcf" > "${OUT}/${CHROM}.popped.vcf"
    cx "$PGGB_SIF" vcfwave -I 1000 -t "$THREADS" "${OUT}/${CHROM}.popped.vcf" | gzip > "${OUT}/${CHROM}.sv.vcf.gz"
    n_raw=$(grep -vc '^#' "${OUT}/${CHROM}.raw.vcf" || true)
    n_pop=$(grep -vc '^#' "${OUT}/${CHROM}.popped.vcf" || true)
    n_wav=$(zcat "${OUT}/${CHROM}.sv.vcf.gz" | grep -vc '^#' || true)
    rm -f "${OUT}/${CHROM}.raw.vcf" "${OUT}/${CHROM}.popped.vcf"
    echo "  sites: ${n_raw} deconstructed, ${n_pop} after vcfbub, ${n_wav} after vcfwave"
    echo ">>> ${CHROM} done: $(date)"
    exit 0
fi

# ----------------------------- summary ----------------------------------------
if [[ "$STAGE" == "summary" ]]; then
    mkdir -p "${SV_DIR}/summary"
    python3 "${REPO}/13_sv_summary.py" --sv-dir "$SV_DIR" --chroms "$CHROM_LIST" \
        --min-sv "$MIN_SV" --min-called-frac "$MIN_CALLED_FRAC" --equal-n "$EQUAL_N" --seed "$SEED" \
        --out "${SV_DIR}/summary"
    exit 0
fi

echo "Usage: bash $0 submit   |   STAGE=chrom (CHROM=...) / STAGE=summary via sbatch" >&2
exit 1
