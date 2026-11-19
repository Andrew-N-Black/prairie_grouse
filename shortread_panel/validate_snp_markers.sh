#!/bin/bash
#SBATCH -A dewoody
#SBATCH -p cpu
#SBATCH -N 1
#SBATCH -n 8
#SBATCH --mem=12G
#SBATCH -t 4:00:00
#SBATCH -J z_snp_validate
#SBATCH -o logs/14i_validate_%x_%A_%a.log
# =============================================================================
# validate_snp_markers.sh -- check the species-diagnostic Z SNPs (step 22)
# in the higher-coverage shotgun birds, as validate_indel_markers.sh does for the indel markers.
#   * each bird's reads are aligned only to +/-300-bp windows around the SNPs
#     (LEPC reference); alignments with > MAX_NM mismatches are dropped, so
#     diverged copies (W gametologs, paralogs) do not count
#   * every bird is genotyped at every SNP and compared with its species'
#     expected allele; a female (one Z) with two alleles is flagged
#   * the windows are also aligned to the rock ptarmigan genome (Z and W)
#
# Usage (from the folder holding these scripts; logs/ must exist):
#   bash validate_snp_markers.sh submit       # prep -> per-bird array -> summary
# Outputs ($OUT):
#   snp.snp_validation.tsv  one row per SNP x pair: birds called and expected in
#                           each species, unexpected, two alleles (female), pass
#   snp.bird_ancestry.tsv   per bird and pair: fraction of reads with each
#                           species' allele across all SNPs (Z ancestry)
#   snp.refcheck_ptarmigan.tsv  best-hit ptarmigan chromosome for each window
# =============================================================================
set -euo pipefail
MODE=${1:-submit}

# ---- edit -------------------------------------------------------------------
BASE=$CLUSTER_SCRATCH/GROUSE/prairie_grouse
OUT=$BASE/14i_validate_snp_markers
SAMPLES=$BASE/14h_validate_indel_markers/samples.tsv     # sample sheet made by 14h
LEPC_REF=$CLUSTER_SCRATCH/GROUSE/grouse_asm/ref/GCF_026119805.1_pur_lepc_1.0_genomic.fa
PTARMIGAN=$CLUSTER_SCRATCH/GROUSE/grouse_asm/ref/GCF_023343835.1_bLagMut1_primary_genomic.fna
PTARMIGAN_MAP=$CLUSTER_SCRATCH/GROUSE/grouse_asm/ref/GCF_023343835.1_bLagMut1_primary_assembly_report.txt
# step-22 SNP table: one row per SNP with ref/alt, allele_GRPC/LEPC/STGR (R/A)
# and the pairs each SNP separates
SNP_TABLE=$CLUSTER_SCRATCH/GROUSE/z_markers/z_markers.final.tsv
COL_CONTIG=chrom
COL_POS=pos            # 1-based
POS_BASE=1
COL_PAIR=panel_confirmed_pairs   # pairs confirmed in the 506-bird panel (dAF >= 0.95);
                                 # use diagnostic_pairs to test all 5,519 primer-designed SNPs
COL_A=allele_a         # only used for a long-format table (one row per SNP x pair)
COL_B=allele_b
COL_KEEP=""            # optional extra filter column
KEEP_VALUES="1,TRUE,True,yes"
MAX_NM=4               # max mismatches per read alignment
CONDA_ENV=zindel
THREADS=${SLURM_NTASKS:-8}
# -----------------------------------------------------------------------------

if [[ -n $CONDA_ENV ]]; then
  module load conda 2>/dev/null || module load anaconda 2>/dev/null || true
  eval "$(conda shell.bash hook)"; set +u; conda activate "$CONDA_ENV"; set -u
fi
for t in minimap2 samtools python3; do command -v $t >/dev/null || { echo "ERROR: $t not found" >&2; exit 1; }; done
HERE=${SLURM_SUBMIT_DIR:-$(dirname "$(readlink -f "$0")")}
PY=$HERE/validate_snp_markers.py; [[ -f $PY ]] || PY=$HERE/shortread_panel/validate_snp_markers.py
[[ -f $PY ]] || { echo "ERROR: $PY not found" >&2; exit 1; }
[[ -s $SAMPLES ]] || { echo "ERROR: $SAMPLES missing (run validate_indel_markers.sh first)" >&2; exit 1; }
mkdir -p "$OUT"/{counts,bam} "$HERE/logs"

case $MODE in
submit)
  [[ $SNP_TABLE == EDIT_ME || ! -s $SNP_TABLE ]] && { echo "ERROR: set SNP_TABLE and its columns" >&2; exit 1; }
  N=$(($(wc -l < "$SAMPLES") - 1))
  j1=$(sbatch --parsable -q standby "$0" prep)
  j2=$(sbatch --parsable -q standby --array=1-$N -d afterok:$j1 "$0" map)
  j3=$(sbatch --parsable -q standby -n 1 --mem=8G -t 1:00:00 -d afterok:$j2 "$0" summarize)
  echo "submitted: prep $j1, map array $j2 (1-$N), summarize $j3"
  ;;

prep)
  cd "$OUT"
  python3 "$PY" regions --table "$SNP_TABLE" --out-prefix snp --col-contig "$COL_CONTIG" --col-pos "$COL_POS" \
    --pos-offset "$POS_BASE" --col-pair "$COL_PAIR" --col-a "$COL_A" --col-b "$COL_B" \
    --col-keep "$COL_KEEP" --keep-values "$KEEP_VALUES"
  awk '{printf "%s:%d-%d\n", $1, $2 + 1, $3}' snp.regions.bed > snp.regions.txt
  samtools faidx "$LEPC_REF" -r snp.regions.txt > snp.windows.fa
  samtools faidx snp.windows.fa
  # windows vs rock ptarmigan (female reference: has Z and W)
  minimap2 -c -x sr -k 11 -w 3 -t "$THREADS" "$PTARMIGAN" snp.windows.fa > snp.ptarmigan.paf 2>/dev/null
  python3 "$PY" refcheck --paf snp.ptarmigan.paf --regions snp.regions.bed --chrmap "$PTARMIGAN_MAP" \
    --label ptarmigan --out snp.refcheck_ptarmigan.tsv
  ;;

map)
  cd "$OUT"
  i=${SLURM_ARRAY_TASK_ID:?run as an array}
  read -r s sp sex r1 r2 < <(sed -n "$((i + 1))p" "$SAMPLES")
  minimap2 -ax sr -t "$THREADS" snp.windows.fa "$r1" "$r2" 2> /dev/null \
    | samtools view -h -F 0x904 -q 20 - \
    | awk -v m=$MAX_NM 'BEGIN{OFS="\t"} /^@/ {print; next}
        { for (k = 12; k <= NF; k++) if ($k ~ /^NM:i:/) { split($k, a, ":"); if (a[3] > m) next }
          print }' \
    | samtools sort -o bam/$s.bam -
  samtools index bam/$s.bam
  samtools mpileup -l snp.sites_window.bed -f snp.windows.fa -q 20 -Q 20 -B -d 1000 \
      --ff UNMAP,SECONDARY,QCFAIL,DUP bam/$s.bam 2> /dev/null \
    | python3 "$PY" count --out counts/$s.counts.tsv
  echo "$s ($sp, $sex): $(tail -n +2 counts/$s.counts.tsv | wc -l) SNP sites with reads"
  ;;

summarize)
  cd "$OUT"
  python3 "$PY" summarize --samples "$SAMPLES" --sites snp.sites.tsv --counts-dir counts --out-prefix snp
  echo; echo "== Z ancestry of the hybrids and any bird with < 90% of its own species' alleles =="
  awk -F'\t' 'NR==1 || $2=="HYB" || ($2!="NA" && index($4,$2)==1 && $6<0.9) || ($2!="NA" && index($4,$2)>1 && $7<0.9)' \
    snp.bird_ancestry.tsv | (column -t 2>/dev/null || cat)
  ;;
*) echo "unknown mode $MODE (submit | prep | map | summarize)" >&2; exit 1 ;;
esac
