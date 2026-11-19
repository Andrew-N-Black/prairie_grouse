#!/bin/bash
#SBATCH -A dewoody
#SBATCH -p cpu
#SBATCH -N 1
#SBATCH -n 16
#SBATCH --mem=64G
#SBATCH -t 4:00:00
#SBATCH -J z_indel_markers
#SBATCH -o logs/14g_z_indel_markers_%j.log
# =============================================================================
# 22b_z_indel_markers.sh -- species-diagnostic Z-linked INDEL markers that can be
# scored by band length on an agarose gel (a single PCR per marker, no
# sequencing), designed from the 40 Z-bearing haplotype assemblies.
#
# Same logic as step 22 (SNP markers), but for length differences:
#   1  extract chr_Z (and chr_Z pieces) from every haplotype that carries it,
#      chr_W from the six females, and the two Z scaffolds of the LEPC reference
#   2  align each Z to the LEPC Z (minimap2 -cx asm5 --cs) and call variants
#      (paftools.js call)
#   3  cluster indels across haplotypes; keep windows where every haplotype of
#      one species differs in net length from every haplotype of another by
#      >= MIN_DIFF bp, with <= MAX_WITHIN bp variation inside each species, and
#      >= 80% of each species' haplotypes called; drop simple-repeat (STR /
#      homopolymer) indels, which are prone to length homoplasy
#   4  Primer3 on the LEPC reference: target = indel window, all other variable
#      sites in any haplotype excluded from primers, products PRODUCT_RANGE bp
#   5  in-silico PCR (seqkit amplicon, exact primer matches) on every Z
#      haplotype (expected band per species), on the W of each female (must give
#      no product) and on the whole LEPC reference genome (exactly one product);
#      each product is also aligned to the female W (minimap2), because Z-W
#      gametologs may amplify despite primer mismatches
#   6  ranked table: band sizes per species, which pairs are resolvable on a gel
#      (>= GEL_MIN_BP and >= 5% of band size), and pass/fail for every check
#
# Usage (from the directory holding both scripts; logs/ must exist):
#   mkdir -p logs
#   sbatch 22b_z_indel_markers.sh                       all steps, normal queue
#   sbatch -q standby 22b_z_indel_markers.sh            all steps, standby (<= 4 h)
#   STEP=pcr,report sbatch 22b_z_indel_markers.sh       chosen steps (comma list):
#        extract | align | candidates | primers | pcr | report
# Output:  $OUT/z_indel.markers.tsv  (pass_all == 1 are the markers to order)
# =============================================================================
set -euo pipefail

# ---- edit -------------------------------------------------------------------
BASE=$CLUSTER_SCRATCH/GROUSE/prairie_grouse
OUT=$BASE/14g_z_indel_markers
# hap_id <tab> species <tab> assembly fasta (final, uncut assemblies; one row per haplotype)
ASSEMBLIES=$BASE/config/z_indel_assemblies.tsv
LEPC_REF=/scratch/gautschi/blackan/GROUSE/grouse_asm/ref/GCF_026119805.1_pur_lepc_1.0_genomic.fa   # soft-masked LEPC reference (pur_lepc_1.0)
Z_SCAFS="NW_026294758.1 NW_026294813.1"            # LEPC Z scaffolds
Z_REGEX='^chr_Z(_p[0-9]+)?$'                        # Z names in the assemblies
W_REGEX='^chr_W(_p[0-9]+)?$'
Z_MIN_BP=20000000   # a haplotype "carries the Z" only with >= this much chr_Z (full Z ~75 Mb;
                    # W-bearing haplotypes keep small chr_Z fragments that must not count)
W_MIN_BP=1         # any chr_W in a haplotype goes into the W check (males have none;
                    # in females most W sequence sits in the Z-bearing haplotype)
MIN_INDEL=10       # smallest indel considered (bp)
MIN_DIFF=15        # min band-length difference between species (bp)
MAX_WITHIN=2       # max band-length range within a species (bp)
PRODUCT_RANGE=150-450
MAX_PER_PAIR=150   # candidates per species pair sent to Primer3 / in-silico PCR
THREADS=${SLURM_NTASKS:-8}
STEP=${STEP:-all}
# -----------------------------------------------------------------------------
[[ $LEPC_REF == EDIT_ME || ! -s $LEPC_REF ]] && { echo "ERROR: set LEPC_REF" >&2; exit 1; }

# Tools: a conda env is the most reliable source (bioconda's minimap2 package
# includes paftools.js and k8). Create once on a login node:
#   conda create -y -n zindel -c conda-forge -c bioconda minimap2 k8 samtools seqkit primer3
CONDA_ENV=zindel        # set to "" to use RCAC modules instead
if [[ -n $CONDA_ENV ]]; then
  module load conda 2>/dev/null || module load anaconda 2>/dev/null || true
  eval "$(conda shell.bash hook)"
  set +u; conda activate "$CONDA_ENV"; set -u
else
  # load one at a time: Lmod aborts the whole list if any one module is missing
  module load biocontainers 2>/dev/null || true
  for m in minimap2 samtools seqkit primer3; do module load $m 2>/dev/null || echo "WARN: no module $m" >&2; done
fi
for t in minimap2 samtools seqkit primer3_core; do
  command -v $t >/dev/null || { echo "ERROR: $t not found" >&2; exit 1; }
done
PAFTOOLS=$(command -v paftools.js || true)
if [[ -z $PAFTOOLS ]]; then
  PAFTOOLS="k8 $(dirname "$(readlink -f "$(command -v minimap2)")")/paftools.js"
fi
echo | $PAFTOOLS version >/dev/null 2>&1 || { echo "ERROR: paftools.js not usable ($PAFTOOLS)" >&2; exit 1; }
echo "Tools: $(minimap2 --version) minimap2; $(samtools --version | head -1); $(seqkit version); paftools: $PAFTOOLS"
# under sbatch $0 is a spooled copy, so look next to where the job was submitted
PY=${SLURM_SUBMIT_DIR:-$(dirname "$(readlink -f "$0")")}/22b_z_indel_markers.py
[[ -f $PY ]] || { echo "ERROR: $PY not found; submit from the directory holding both scripts" >&2; exit 1; }
mkdir -p "$OUT"/{z,w,calls,pcr,logs} logs
cd "$OUT"
run() { [[ $STEP == all || ",$STEP," == *",$1,"* ]]; }   # STEP may be a comma list

# ---- 1 extract --------------------------------------------------------------
if run extract; then
  samtools faidx "$LEPC_REF"
  samtools faidx "$LEPC_REF" $Z_SCAFS > ref_Z.fa            # keeps soft-masking
  samtools faidx ref_Z.fa
  : > sex_chrom_lengths.tsv
  rm -f z/*.fa w/*.fa
  printf "hap_id\tspecies\tz_fasta\tcalls\n" > haplotypes.tsv
  while IFS=$'\t' read -r hap sp fa; do
    [[ $hap == hap_id || -z $hap ]] && continue
    [[ -f $fa.fai ]] || samtools faidx "$fa"
    zlen=$(awk -v r="$Z_REGEX" '$1 ~ r {s+=$2} END {print s+0}' "$fa.fai")
    wlen=$(awk -v r="$W_REGEX" '$1 ~ r {s+=$2} END {print s+0}' "$fa.fai")
    z=$(awk -v r="$Z_REGEX" '$1 ~ r {print $1}' "$fa.fai")
    w=$(awk -v r="$W_REGEX" '$1 ~ r {print $1}' "$fa.fai")
    printf "%s\t%s\tchr_Z=%d\tchr_W=%d\n" "$hap" "$sp" "$zlen" "$wlen" >> sex_chrom_lengths.tsv
    if (( wlen >= W_MIN_BP )); then samtools faidx "$fa" $w > w/$hap.fa; fi
    if (( zlen >= Z_MIN_BP )); then
      samtools faidx "$fa" $z > z/$hap.fa
      printf "%s\t%s\t%s\t%s\n" "$hap" "$sp" "$OUT/z/$hap.fa" "$OUT/calls/$hap.txt" >> haplotypes.tsv
    fi
  done < "$ASSEMBLIES"
  echo "Z-bearing haplotypes: $(($(wc -l < haplotypes.tsv) - 1)); haplotypes with chr_W: $(ls w | wc -l)"
fi

# ---- 2 align + call ---------------------------------------------------------
if run align; then
  tail -n +2 haplotypes.tsv | cut -f1 | while read -r hap; do
    minimap2 -cx asm5 --cs -t "$THREADS" ref_Z.fa z/$hap.fa 2> logs/minimap2_$hap.log \
      | sort -k6,6 -k8,8n \
      | $PAFTOOLS call -s "$hap" - > calls/$hap.txt 2> logs/paftools_$hap.log
  done
fi

# ---- 3 candidate indels -----------------------------------------------------
if run candidates; then
  python3 "$PY" candidates --haps haplotypes.tsv --ref-z ref_Z.fa --out-prefix z_indel \
    --min-indel $MIN_INDEL --min-diff $MIN_DIFF --max-within $MAX_WITHIN
fi

# ---- 4 primers --------------------------------------------------------------
if run primers; then
  python3 "$PY" primers --candidates z_indel.candidates.tsv --ref-z ref_Z.fa \
    --varsites z_indel.all_variable_sites.tsv --out-prefix z_indel \
    --product-range $PRODUCT_RANGE --max-per-pair $MAX_PER_PAIR
fi

# ---- 5 in-silico PCR --------------------------------------------------------
if run pcr; then
  P=z_indel.primers_seqkit.tsv
  pcr() { seqkit amplicon -p "$P" --bed -j 2 --quiet "$2" | cut -f1-6 > "pcr/$1.bed"; }
  export -f pcr; export P
  { for f in z/*.fa; do h=$(basename $f .fa); echo "Z.$h $f"; done
    for f in w/*.fa; do [[ -e $f ]] || continue; h=$(basename $f .fa); echo "W.$h $f"; done
  } | xargs -P $((THREADS / 2)) -n 2 bash -c 'pcr "$0" "$1"'
  seqkit amplicon -p "$P" --bed -j "$THREADS" --quiet "$LEPC_REF" | cut -f1-6 > pcr/refgenome.bed
  # exact primer matches miss diverged Z-W gametologs, so also align each LEPC
  # product to the W of every female (>= 80% identity over >= 50% = W homolog)
  if compgen -G "w/*.fa" > /dev/null; then
    cat w/*.fa > pcr/W_all.fa
    minimap2 -c -x sr -k 11 -w 3 -N 50 -t "$THREADS" pcr/W_all.fa z_indel.amplicons.fa \
      > pcr/W_homology.paf 2> logs/minimap2_W.log
  fi
fi

# ---- 6 report ---------------------------------------------------------------
if run report; then
  python3 "$PY" report --primers z_indel.primers.tsv --haps haplotypes.tsv \
    --pcr-dir pcr --out-prefix z_indel
  echo "Markers to order: awk -F'\t' 'NR==1 || \$6==1' $OUT/z_indel.markers.tsv"
fi
