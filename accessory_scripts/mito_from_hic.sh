#!/bin/bash
# =============================================================================
# accessory_scripts/mito_from_hic.sh — mitogenome from the Hi-C (Illumina) reads,
# for birds where MitoHiFi cannot assemble one from HiFi (STGR F5485, F5502).
#
# WHY: in these two birds the HiFi reads that match the mitogenome appear to be
# mostly NUMT-bearing nuclear reads. Both MitoHiFi runs built mostly nuclear
# contigs that cover only 8-29% of the mitogenome, whatever the read-length
# cutoff. Short-read libraries such as Hi-C carry mtDNA at hundreds to thousands
# of times nuclear depth, so a consensus called from them is dominated by the
# true mitogenome rather than NUMTs.
#
# HOW: reads are mapped to a PASS mitogenome of the same species (STGR F5457 by
# default); haploid variants are called and applied to give the bird's own
# consensus; sites under MIN_DEPTH are masked to N; MitoHiFi then re-annotates
# and rotates the result (-c).
# CAVEAT: a reference-guided consensus can under-call large indels or
# rearrangements, mainly in the control region. Report these two as
# "Hi-C read consensus" in the methods; the other 21 are de novo HiFi assemblies.
#
# RUN (rows of assembly_manifest.tsv; F5485 = 7, F5502 = 8):
#   sbatch -q standby -t 04:00:00 --array=7,8 accessory_scripts/mito_from_hic.sh
# VALIDATE the method on a bird whose HiFi mitogenome passed (F5457 = row 0), using
# a different guide; the result goes to STGR_F5457.hic_test and is compared afterwards:
#   sbatch -q standby -t 04:00:00 --array=0 --export=ALL,OUT_SUFFIX=.hic_test,GUIDE=$ASM/mitogenomes/STGR_F5462/final_mitogenome.fasta accessory_scripts/mito_from_hic.sh
# Then: bash accessory_scripts/mito_qc.sh qc
# =============================================================================
#SBATCH -A dewoody
#SBATCH -p cpu
#SBATCH -J mito_hic
#SBATCH -N 1
#SBATCH -n 1
#SBATCH -c 16
#SBATCH --mem=32G
#SBATCH -t 04:00:00
#SBATCH -o logs/%x_%A_%a.out
#SBATCH -e logs/%x_%A_%a.err

set -euo pipefail
unset LD_PRELOAD
module load biocontainers minimap2 samtools bcftools

PROJECT_DIR="${CLUSTER_SCRATCH}/GROUSE/grouse_asm"
MANIFEST="${SLURM_SUBMIT_DIR:-$PWD}/assembly_manifest.tsv"
MITO_DIR="${PROJECT_DIR}/mitogenomes"
REF_DIR="${MITO_DIR}/reference"
SIF="${PROJECT_DIR}/containers/mitohifi.sif"
GUIDE="${GUIDE:-${MITO_DIR}/STGR_F5457/final_mitogenome.fasta}"   # same-species PASS mitogenome
MIN_DEPTH="${MIN_DEPTH:-10}"
THREADS="${SLURM_CPUS_PER_TASK:-8}"
APPTAINER_BIN="$(command -v apptainer || command -v singularity)"

LINE=$(grep -v '^#' "$MANIFEST" | tail -n +2 | sed -n "$((SLURM_ARRAY_TASK_ID + 1))p")
IFS=$'\t' read -r SAMPLE SPECIES _HIFI HIC_R1 HIC_R2 _ONT <<< "$LINE"
SPECIES="${SPECIES//[[:space:]]/}"
OUT="${MITO_DIR}/${SPECIES}_${SAMPLE}${OUT_SUFFIX:-}"   # OUT_SUFFIX=.hic_test leaves an existing result untouched

# keep any earlier attempt
if [[ -d "$OUT" && ! -f "$OUT/HIC_CONSENSUS" ]]; then
    mv "$OUT" "${OUT}.failed_hifi_$(date +%s)"
fi
mkdir -p "$OUT/hic"; cd "$OUT/hic"
echo ">>> ${SPECIES} ${SAMPLE}: Hi-C consensus against $(basename "$(dirname "$GUIDE")")"

sed "1s/.*/>guide/" "$GUIDE" > guide.fa
samtools faidx guide.fa
IFS=',' read -ra R1S <<< "$HIC_R1"
IFS=',' read -ra R2S <<< "$HIC_R2"

# 1. map Hi-C reads (primary, MAPQ>=20). Hi-C pairs can be chimeric across
#    loci, so properly-paired is not required; each read stands on its own.
minimap2 -t "$THREADS" -ax sr guide.fa <(zcat -f "${R1S[@]}") <(zcat -f "${R2S[@]}") 2>/dev/null \
  | samtools view -b -F 0x904 -q 20 -e '[NM] <= 3' - | samtools sort -@ 4 -o hic_mito.bam -
# [NM] <= 3: a true mitochondrial read differs from a same-species guide at well
# under 1% of sites; NUMT copies have diverged several percent, so most NUMT reads
# carry more than 3 mismatches per ~150 bp and are dropped here.
samtools index hic_mito.bam
samtools depth -a hic_mito.bam > depth.txt
awk '{s+=$3; if($3<m||NR==1)m=$3} END{printf ">>> depth: mean %.0fx, min %dx\n", s/NR, m}' depth.txt

# 2. haploid calls -> consensus; mask low-depth sites
bcftools mpileup -f guide.fa -d 100000 -Q 20 -q 20 -a AD,DP -Ob -o pileup.bcf hic_mito.bam
bcftools call --ploidy 1 -mv -Oz -o calls.vcf.gz pileup.bcf
# mixed sites: depth >= 10 but the commonest base is < 80% of reads. A handful is
# noise; dozens mean NUMT reads are still in the mix and the consensus is unsafe.
MIXED=$(bcftools query -f '[%AD]\n' pileup.bcf | awk -F',' '{s=0;m=0; for(i=1;i<=NF;i++){s+=$i; if($i>m)m=$i}} s>=10 && m/s<0.8 {n++} END{print n+0}')
echo ">>> mixed sites (major allele < 80%, depth >= 10): ${MIXED}"
bcftools index calls.vcf.gz
awk -v m="$MIN_DEPTH" '$3<m {print $1"\t"$2-1"\t"$2}' depth.txt > lowdepth.bed
bcftools consensus -f guide.fa -m lowdepth.bed calls.vcf.gz | sed "1s/.*/>${SPECIES}_${SAMPLE}_mito_hic/" > consensus.fa
echo ">>> variants vs guide: $(bcftools view -H calls.vcf.gz | wc -l)   masked (<${MIN_DEPTH}x): $(wc -l < lowdepth.bed) bp"

# 3. annotate + rotate with MitoHiFi (contig mode), promote to the sample folder
cd "$OUT"
"$APPTAINER_BIN" exec --cleanenv --bind "${CLUSTER_SCRATCH}" "$SIF" \
  mitohifi.py -c hic/consensus.fa -f "${REF_DIR}/ref.fasta" -g "${REF_DIR}/ref.gb" -t 4 -o 2 -a animal
cp hic/hic_mito.bam* . 2>/dev/null || true
touch HIC_CONSENSUS
echo ">>> final length: $(grep -v '>' final_mitogenome.fasta | tr -d '\n' | wc -c) bp"
